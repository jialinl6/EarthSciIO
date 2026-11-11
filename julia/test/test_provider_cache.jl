# The Provider's last-sample memo (`materialize`: "the result is shared"). It may
# change how often a blob is READ and nothing else, so every sample below is
# compared, field for field with `isequal` on every array, against the same
# provider asked with the memo off — the path every sample took before it.
# Reuses `same_native` from test_provider.jl (included first by runtests.jl).

# A netcdf reader that counts its decodes, which is what the memo removes. It
# declares the real reader's options so the record pushdown still applies.
struct CountingNetCDFReader <: EarthSciIO.Reader end
const NETCDF_DECODES = Threads.Atomic{Int}(0)
function EarthSciIO.read_native(::CountingNetCDFReader, path::AbstractString;
                                variables = nothing, select = nothing, records = nothing)
    Threads.atomic_add!(NETCDF_DECODES, 1)
    return read_native(FORMAT_REGISTRY["netcdf"], path; variables = variables,
                       select = select, records = records)
end
EarthSciIO.supports_selection(::CountingNetCDFReader) = true
register!(FORMAT_REGISTRY, "netcdf-counting", CountingNetCDFReader())
decodes() = NETCDF_DECODES[]

memo_off(f) = withenv(f, "EARTHSCIIO_PROVIDER_CACHE" => "0")

# Same names, dims, attrs and eltypes, and `isequal` on every array (NaN == NaN),
# on top of `same_native`'s bit-pattern equality.
function identical_native(a, b)
    same_native(a, b) || return false
    for n in vcat(variable_names(a), coord_names(a))
        isequal(a[n].data, b[n].data) || return false
    end
    return true
end

# Pin a file's mtime without touching its bytes. A just-written file is "racy"
# (its mtime is inside the 2 s margin), and no cache keeps anything read from it.
function set_mtime!(path, t)
    fh = Base.Filesystem.open(path, Base.Filesystem.JL_O_WRONLY)
    try
        Base.Filesystem.futime(fh, Float64(t), Float64(t))
    finally
        close(fh)
    end
    return path
end
backdate!(path, seconds_ago = 60) = set_mtime!(path, floor(time()) - seconds_ago)
# Inside the margin for the whole test, however long it runs: ahead of the clock.
make_racy!(path) = set_mtime!(path, floor(time()) + 3600)

# A daily file of `nrec` hourly records on a 3 x 4 grid. `t2m` carries the record
# axis, varies over space (so a `select` window shows up as different numbers)
# and has one masked cell; `orog` has no record axis and must pass through.
function write_hourly_nc(path; nrec::Int, base::Float64, t0::Float64)
    NC = EarthSciIO.NCDatasets
    NC.NCDataset(path, "c") do ds
        NC.defDim(ds, "time", nrec)
        NC.defDim(ds, "latitude", 3)
        NC.defDim(ds, "longitude", 4)
        tv = NC.defVar(ds, "time", Float64, ("time",))
        tv.attrib["units"] = "hours since 2018-11-09 00:00:00"
        tv.attrib["calendar"] = "standard"
        tv[:] = collect(t0:(t0 + nrec - 1))
        NC.defVar(ds, "latitude", Float64, ("latitude",))[:] = [40.0, 41.0, 42.0]
        NC.defVar(ds, "longitude", Float64, ("longitude",))[:] = [-90.0, -89.0, -88.0, -87.0]
        v = NC.defVar(ds, "t2m", Float64, ("longitude", "latitude", "time");
                      fillvalue = -9999.0)
        for r in 1:nrec, la in 1:3, lo in 1:4
            v[lo, la, r] = base + r + la / 10 + lo / 100
        end
        v[2, 3, 1] = missing
        NC.defVar(ds, "orog", Float64, ("longitude", "latitude"))[:, :] .= 7.0
    end
    return path
end

blob_path(root, url) = joinpath(root, "v1", "blobs", cache_key(url)[1:2], cache_key(url) * ".nc")

# Stage `url => (nrec, base, t0)` blobs into a fresh offline cache, safely old.
function stage_hourly(specs)
    root = mktempdir()
    for (u, (nrec, base, t0)) in specs
        p = blob_path(root, u)
        mkpath(dirname(p))
        backdate!(write_hourly_nc(p; nrec = nrec, base = base, t0 = t0))
    end
    return Cache(LocalStore(root); offline = true, verify = false), root
end

# Three daily files of four hourly records: a 12-tick cadence whose brackets
# cross a file seam at ticks 4 and 8.
const DAY_URLS = ["https://data.earthsci.dev/memo/day1.nc",
                  "https://data.earthsci.dev/memo/day2.nc",
                  "https://data.earthsci.dev/memo/day3.nc"]
day_url(t) = DAY_URLS[clamp(Int(fld(t, 4.0)) + 1, 1, 3)]
const DAY_TICKS = collect(0.0:11.0)
stage_days() = stage_hourly([DAY_URLS[d] => (4, 100.0 * d, 4.0 * (d - 1)) for d in 1:3])

@testset "Provider last-sample memo" begin
    c, root = stage_days()
    WIN = Dict("axes" => Any["all", Dict("indices" => [2, 0]), Dict("slice" => [1, 3])])

    # Sample `p` at each of `probes` with the memo on and check the memo's whole
    # contract: a repeated tick reads nothing, a moved tick reads, and either way
    # the sample is identical to the memo-off read. (The first probe may hit what
    # an earlier sweep left, so its read count is not asserted.)
    function sweep(p, probes; select = nothing)
        prev = nothing
        for t in probes
            before = decodes()
            got = materialize(p, t; select = select)
            tick = EarthSciIO._tick_index(p, t)
            if prev !== nothing
                tick == prev ? (@test decodes() == before) : (@test decodes() > before)
            end
            prev = tick
            @test identical_native(got, memo_off(() -> materialize(p, t; select = select)))
        end
    end
    forward  = [0.0, 0.5, 1.0, 1.0, 2.0, 3.0, 3.5, 4.0, 5.0, 7.0, 7.5, 8.0, 11.0, 11.0, 30.0]
    backward = [11.0, 10.5, 10.0, 8.0, 7.5, 7.0, 4.0, 3.5, 3.0, 1.0, 0.0, 0.0]
    repeated = [5.0, 5.0, 5.25, 5.0, 5.999]

    @testset "memo on == memo off: rps=$(repr(rps)), pushdown=$push, select=$(sel !== nothing)" for
            rps in (nothing, 1, 2), push in (true, false), sel in (nothing, WIN)
        kw = push ? NamedTuple() : (reader_kwargs = (records = nothing,),)
        p = discrete_provider(c, day_url, DAY_TICKS; format = "netcdf-counting",
                              time_dim = "time", records_per_sample = rps, kw...)
        @test EarthSciIO._memo_applies(p)
        @test EarthSciIO._records_pushdown(p) == push
        sweep(p, forward; select = sel)
        sweep(p, backward; select = sel)        # the adjoint replays the forcing backwards
        sweep(p, repeated; select = sel)
    end

    @testset "a baked select is the key's select too" begin
        p = discrete_provider(c, day_url, DAY_TICKS; format = "netcdf-counting",
                              time_dim = "time", records_per_sample = 2,
                              reader_kwargs = (select = WIN,))
        sweep(p, vcat(forward, backward))
        # A per-call select overrides the baked one, and is a different key.
        a = materialize(p, 2.0)
        n = decodes()
        b = materialize(p, 2.0; select = Dict("axes" => Any["all", "all", "all"]))
        @test decodes() > n
        @test size(a["t2m"].data) == (2, 2, 2) && size(b["t2m"].data) == (2, 3, 4)
    end

    @testset "the numbers, so on == off cannot pass with both wrong" begin
        p = discrete_provider(c, day_url, DAY_TICKS; format = "netcdf",
                              time_dim = "time", records_per_sample = 2)
        cell(t) = materialize(p, t)["t2m"].data[:, 1, 1]   # lat 40, lon -90
        @test cell(1.0) ≈ [102.11, 103.11]      # interior
        @test cell(3.0) ≈ [104.11, 201.11]      # the day1 -> day2 seam
        @test cell(3.0) ≈ [104.11, 201.11]      # ...again, from the memo
        @test cell(0.0) ≈ [101.11, 102.11]      # backwards
        @test cell(11.0) ≈ [304.11, 304.11]     # end of data degenerates
        @test isnan(materialize(p, 0.0)["t2m"].data[1, 3, 2])    # the masked cell
    end

    @testset "same tick returns the SAME object; the memo holds one sample" begin
        p = discrete_provider(c, day_url, DAY_TICKS; format = "netcdf",
                              time_dim = "time", records_per_sample = 2)
        a = materialize(p, 5.0)
        @test materialize(p, 5.5) === a          # one tick, any t inside it
        @test refresh(p, 5.0) === a              # refresh is materialize
        b = materialize(p, 6.0)
        @test b !== a
        @test p.memo.value === b && p.memo.key[1] == 7   # only the last sample is kept
        @test materialize(p, 5.0) !== a          # going back re-reads; nothing older is kept
    end

    @testset "select: equal content hits, a different window misses" begin
        p = discrete_provider(c, day_url, DAY_TICKS; format = "netcdf-counting",
                              time_dim = "time", records_per_sample = 2)
        sel = deepcopy(WIN)
        materialize(p, 1.0; select = sel)
        n = decodes()
        @test identical_native(materialize(p, 1.0; select = deepcopy(WIN)),
                               memo_off(() -> materialize(p, 1.0; select = WIN)))
        @test decodes() == n + 1                 # the hit read nothing; memo_off read once
        # Editing the caller's Dict in place is a new key, not a stale hit: the memo
        # keeps a copy of the select it was read under.
        sel["axes"][2] = Dict("indices" => [1])
        n = decodes()
        w = materialize(p, 1.0; select = sel)
        @test decodes() > n
        @test size(w["t2m"].data) == (2, 1, 2)
        @test identical_native(w, memo_off(() -> materialize(p, 1.0; select = sel)))
    end

    @testset "a re-fetched, rewritten or touched blob is re-read" begin
        c2, root2 = stage_days()
        p = discrete_provider(c2, day_url, DAY_TICKS; format = "netcdf-counting",
                              time_dim = "time", records_per_sample = 2)
        path = blob_path(root2, DAY_URLS[1])
        @test materialize(p, 1.0)["t2m"].data[:, 1, 1] ≈ [102.11, 103.11]
        n = decodes()
        materialize(p, 1.0)
        @test decodes() == n                     # warm

        # touched: same bytes, new mtime -> re-read, same numbers
        backdate!(path, 30)
        @test materialize(p, 1.0)["t2m"].data[:, 1, 1] ≈ [102.11, 103.11]
        @test decodes() == n + 1

        # replaced by a rename (how the cache commits a re-fetched blob), same size
        tmp = write_hourly_nc(tempname() * ".nc"; nrec = 4, base = 500.0, t0 = 0.0)
        @test filesize(tmp) == filesize(path)    # size alone cannot betray it
        mv(tmp, path; force = true)
        backdate!(path, 45)
        @test materialize(p, 1.0)["t2m"].data[:, 1, 1] ≈ [502.11, 503.11]
        @test decodes() == n + 2

        # rewritten in place (same inode, same size) with an older-looking mtime
        tmp = write_hourly_nc(tempname() * ".nc"; nrec = 4, base = 700.0, t0 = 0.0)
        write(path, read(tmp))
        backdate!(path, 90)
        @test materialize(p, 1.0)["t2m"].data[:, 1, 1] ≈ [702.11, 703.11]
        @test decodes() == n + 3
        materialize(p, 1.0)
        @test decodes() == n + 3                 # and kept again once it is safely old
    end

    @testset "a blob modified within the racy margin is never kept" begin
        c3, root3 = stage_days()
        path = blob_path(root3, DAY_URLS[2])
        make_racy!(path)
        p = discrete_provider(c3, day_url, DAY_TICKS; format = "netcdf-counting",
                              time_dim = "time", records_per_sample = 2)
        n = decodes()
        a = materialize(p, 5.0)
        b = materialize(p, 5.0)
        @test decodes() == n + 2                 # read both times
        @test p.memo.value === nothing
        @test identical_native(a, b)
        # The bracket at tick 4 also reads day2 (its successor), so it is not kept either.
        materialize(p, 3.0); materialize(p, 3.0)
        @test decodes() == n + 6
        backdate!(path)
        materialize(p, 5.0)
        m = decodes()
        @test materialize(p, 5.0) === p.memo.value
        @test decodes() == m
    end

    @testset "EARTHSCIIO_PROVIDER_CACHE=0 restores the read-every-time path" begin
        p = discrete_provider(c, day_url, DAY_TICKS; format = "netcdf-counting",
                              time_dim = "time", records_per_sample = 2)
        for v in ("0", "false", "no", "off", " OFF ")
            withenv("EARTHSCIIO_PROVIDER_CACHE" => v) do
                @test !EarthSciIO._memo_applies(p)
                n = decodes()
                a = materialize(p, 2.0)
                b = materialize(p, 2.0)
                @test decodes() == n + 2
                @test a !== b && identical_native(a, b)
            end
        end
        for v in ("1", "true", "", "yes")
            withenv(() -> (@test EarthSciIO._memo_applies(p)),
                    "EARTHSCIIO_PROVIDER_CACHE" => v)
        end
    end

    @testset "only time_dim providers over a whole-file reader" begin
        era5 = "https://data.earthsci.dev/era5/2018/11/20181108.nc"
        corpus = Cache(LocalStore(joinpath(CORPUS, "cache")); offline = true, verify = true)
        @test !EarthSciIO._memo_applies(const_provider(corpus, era5; format = "netcdf"))
        @test !EarthSciIO._memo_applies(discrete_provider(corpus, era5, [0.0, 1.0];
                                                          format = "netcdf"))
        @test !EarthSciIO._memo_applies(discrete_provider(corpus, era5, [0.0, 1.0];
                                                          format = "zarr", time_dim = "time"))
        @test EarthSciIO._memo_applies(discrete_provider(corpus, era5, [0.0, 1.0];
                                                         format = "netcdf", time_dim = "time"))
    end

    @testset "concurrent tasks sampling one provider get their own tick's sample" begin
        p = discrete_provider(c, day_url, DAY_TICKS; format = "netcdf",
                              time_dim = "time", records_per_sample = 2)
        probes = repeat([0.0, 3.0, 3.5, 4.0, 11.0, 7.0], 8)
        want = Dict(t => memo_off(() -> materialize(p, t)) for t in unique(probes))
        got = Vector{Any}(undef, length(probes))
        @sync for (i, t) in enumerate(probes)
            Threads.@spawn got[i] = materialize(p, t)
        end
        @test all(identical_native(got[i], want[t]) for (i, t) in enumerate(probes))
    end
end
