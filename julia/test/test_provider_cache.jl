# The Provider's last-sample memo (`materialize`: "the result is shared"). It may
# change how often a blob is READ and nothing else, so every sample below is
# compared, field for field with `isequal` on every array, against the same
# provider asked with the memo off — the path every sample took before it.
# Reuses `same_native` from test_provider.jl (included first by runtests.jl).
using Dates

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

# The same as a `temporal` provider (spec/esm-sources.md §2): three daily files of
# 24 hourly records from 2018-11-09, located on the calendar and by each file's
# own time axis; the URL is a function of the file anchor in Unix seconds. Its
# times are seconds from `origin` (the Unix epoch when `nothing`).
const T0 = DateTime(2018, 11, 9)
const TDAY_URLS = ["https://data.earthsci.dev/memo/temporal/day$d.nc" for d in 1:3]
tday_url(u) = TDAY_URLS[clamp(Int(fld(u - datetime2unix(T0), 86400.0)) + 1, 1, 3)]
stage_tdays(days = 1:3) =
    stage_hourly([TDAY_URLS[d] => (24, 100.0 * d, 24.0 * (d - 1)) for d in days])
hours(h; origin = T0) = (origin === nothing ? datetime2unix(T0) : 0.0) + 3600.0 * h
function temporal_provider(
        c; rps = 2, origin = T0, window_h = (0, 72), format = "netcdf-counting", kw...
    )
    s = SourceTemporal(
        T0, "PT1H", "P1D"; stop = T0 + Day(3), time_dim = "time", records_per_sample = rps
    )
    return Provider(
        c, tday_url; format = format, temporal = s, time_origin = origin,
        window = (hours(window_h[1]; origin = origin), hours(window_h[2]; origin = origin)), kw...
    )
end

@testset "Provider last-sample memo" begin
    c, root = stage_days()
    WIN = Dict("axes" => Any["all", Dict("indices" => [2, 0]), Dict("slice" => [1, 3])])

    # Sample `p` at each of `probes` with the memo on and check the memo's whole
    # contract: a repeated record reads nothing, a moved record reads, and either
    # way the sample is identical to the memo-off read. (The first probe may hit
    # what an earlier sweep left, so its read count is not asserted.) The record
    # is the cadence tick, or a `temporal` provider's record anchor.
    function sweep(p, probes; select = nothing)
        prev = nothing
        for t in probes
            before = decodes()
            got = materialize(p, t; select = select)
            tick = EarthSciIO._memo_record(p, t)
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
        # A `temporal` provider always has a time_dim, and a windowed one too.
        tp = temporal_provider(c; format = "netcdf")
        @test EarthSciIO._memo_applies(tp)
        @test EarthSciIO._memo_applies(with_window(tp, DimWindow("latitude" => (0, 2))))
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

@testset "Provider last-sample memo — temporal providers (EarthSciIO 0.4)" begin
    c, root = stage_tdays()
    function sweep(p, probes; select = nothing)
        prev = nothing
        for t in probes
            before = decodes()
            got = materialize(p, t; select = select)
            rec = EarthSciIO._memo_record(p, t)
            if prev !== nothing
                rec == prev ? (@test decodes() == before) : (@test decodes() > before)
            end
            prev = rec
            @test identical_native(got, memo_off(() -> materialize(p, t; select = select)))
        end
    end
    # Hours since 2018-11-09: records inside a day, both day seams, the last record
    # (whose bracket degenerates: its successor is at `stop`), and back again.
    fwd = [0.0, 0.5, 1.0, 1.0, 1.75, 22.0, 23.0, 23.5, 24.0, 25.0, 47.0, 47.99, 48.0, 71.0, 71.5]
    bwd = reverse(fwd)
    rep = [30.0, 30.25, 30.0, 30.999, 30.5]

    @testset "memo on == memo off: rps=$rps, pushdown=$push, origin=$(origin === nothing ? "unix" : "T0")" for
            rps in (1, 2), push in (true, false), origin in (T0, nothing)
        kw = push ? NamedTuple() : (reader_kwargs = (records = nothing,),)
        p = temporal_provider(c; rps = rps, origin = origin, kw...)
        @test EarthSciIO._memo_applies(p)
        @test EarthSciIO._records_pushdown(p) == push
        for probes in (fwd, bwd, rep)
            sweep(p, [hours(h; origin = origin) for h in probes])
        end
    end

    @testset "the numbers, the time coordinate and the origin" begin
        p = temporal_provider(c; format = "netcdf")
        cell(h) = materialize(p, hours(h))["t2m"].data[:, 1, 1]   # lat 40, lon -90
        @test cell(1.0) ≈ [102.11, 103.11]              # interior: records 2 and 3 of day 1
        @test cell(23.0) ≈ [124.11, 201.11]             # the day1 -> day2 seam
        @test cell(23.5) ≈ [124.11, 201.11]             # ...again, from the memo
        @test cell(71.0) ≈ [324.11, 324.11]             # last record: held
        @test isnan(materialize(p, hours(0))["t2m"].data[1, 3, 2])  # the masked cell
        # The bracket's time coordinate counts from the provider's origin, which
        # is in the key with the cadence: a provider on the Unix clock returns the
        # same data at the same instant, with its own coordinate and units.
        u = temporal_provider(c; format = "netcdf", origin = nothing)
        for h in (0.0, 5.5, 23.0, 47.0)
            a, b = materialize(p, hours(h)), materialize(u, hours(h; origin = nothing))
            @test isequal(a["t2m"].data, b["t2m"].data)
            @test a.coords["time"].data == [3600.0 * floor(h), 3600.0 * (floor(h) + 1)]
            @test b.coords["time"].data == a.coords["time"].data .+ datetime2unix(T0)
            @test a.coords["time"].attrs["units"] == "seconds since 2018-11-09T00:00:00Z"
            @test b.coords["time"].attrs["units"] == "seconds since 1970-01-01T00:00:00Z"
        end
        @test p.memo.key[1] == T0 + Hour(47)           # the record anchor
        @test p.memo.key[2] == (p.plan.temporal, p.plan.clock)
        @test u.memo.key[2][2].origin === nothing && p.memo.key[2][2].origin == T0
    end

    @testset "the run window chooses refresh times, never a sample" begin
        p = temporal_provider(c)
        q = temporal_provider(c; window_h = (24, 48))
        @test refresh_times(p) != refresh_times(q)
        for h in (0.0, 23.0, 30.0, 47.0, 60.0)          # inside q's window and outside it
            @test identical_native(materialize(q, hours(h)), materialize(p, hours(h)))
            @test identical_native(materialize(q, hours(h)), memo_off(() -> materialize(q, hours(h))))
        end
    end

    @testset "a successor that cannot be read holds the bracket, and is not kept" begin
        c2, _ = stage_tdays(1:2)                        # day 3 is not in the offline cache
        p = temporal_provider(c2)
        want = memo_off(() -> materialize(p, hours(47)))
        @test want["t2m"].data[1, 1, 1] == want["t2m"].data[2, 1, 1]       # held
        @test want.coords["time"].data[1] == want.coords["time"].data[2]
        n = decodes()
        a = materialize(p, hours(47))
        b = materialize(p, hours(47.5))
        @test identical_native(a, want) && identical_native(b, want)
        @test decodes() == n + 2                        # read both times
        @test p.memo.value === nothing
        materialize(p, hours(46))                       # a readable bracket is kept again
        m = decodes()
        @test materialize(p, hours(46.5)) === p.memo.value && decodes() == m
    end

    @testset "a re-fetched, rewritten or touched blob is re-read" begin
        c3, root3 = stage_tdays()
        p = temporal_provider(c3)
        path = blob_path(root3, TDAY_URLS[2])
        cell() = materialize(p, hours(30))["t2m"].data[:, 1, 1]
        @test cell() ≈ [207.11, 208.11]
        n = decodes()
        cell()
        @test decodes() == n
        backdate!(path, 30)
        @test cell() ≈ [207.11, 208.11] && decodes() == n + 1
        mv(write_hourly_nc(tempname() * ".nc"; nrec = 24, base = 500.0, t0 = 24.0), path;
           force = true)
        backdate!(path, 45)
        @test cell() ≈ [507.11, 508.11] && decodes() == n + 2
        make_racy!(path)
        cell(); cell()
        @test decodes() == n + 4 && p.memo.value === nothing
    end

    @testset "named-dimension windows: with_window and a per-call DimWindow" begin
        p = temporal_provider(c)
        w = DimWindow("latitude" => (0, 2), "longitude" => [3, 1])
        q = with_window(p, w)
        @test q.memo !== p.memo                         # a windowed provider keeps its own
        sweep(q, [hours(h) for h in vcat(fwd, bwd)])
        @test size(materialize(q, hours(5))["t2m"].data) == (2, 2, 2)
        # Per call: an equal window built afresh is the same key, a different one is not.
        materialize(p, hours(5); select = w)
        n = decodes()
        r = materialize(p, hours(5); select = DimWindow("latitude" => (0, 2), "longitude" => [3, 1]))
        @test decodes() == n
        @test identical_native(r, memo_off(() -> materialize(q, hours(5))))
        s = materialize(p, hours(5); select = DimWindow("latitude" => (1, 3)))
        @test decodes() > n
        @test identical_native(s, memo_off(() -> materialize(p, hours(5); select = DimWindow("latitude" => (1, 3)))))
        # sample_windows_at still equals each window sampled on its own.
        ws = [w, DimWindow("latitude" => [2])]
        got = sample_windows_at(p, hours(7), ws)
        for (i, wi) in enumerate(ws)
            @test identical_native(got[i], memo_off(() -> sample_at(with_window(p, wi), hours(7))))
        end
    end

    @testset "a PrefetchGroup serves what it serves without the memo" begin
        function epochs(memo)
            return withenv("EARTHSCIIO_PROVIDER_CACHE" => memo) do
                g = PrefetchGroup(Dict("m" => temporal_provider(c; window_h = (0, 30))); depth = 2)
                try
                    start!(g, :forward)
                    [take!(g, k)["m"].data for k in eachindex(boundaries(g))]
                finally
                    close(g)
                end
            end
        end
        on, off = epochs("1"), epochs("0")
        @test length(on) == length(off) >= 30
        @test all(identical_native(a, b) for (a, b) in zip(on, off))
    end

    @testset "concurrent tasks sampling one temporal provider" begin
        p = temporal_provider(c; format = "netcdf")
        probes = repeat([hours(h) for h in (0.0, 23.0, 23.5, 24.0, 71.0, 47.0)], 8)
        want = Dict(t => memo_off(() -> materialize(p, t)) for t in unique(probes))
        got = Vector{Any}(undef, length(probes))
        @sync for (i, t) in enumerate(probes)
            Threads.@spawn got[i] = materialize(p, t)
        end
        @test all(identical_native(got[i], want[t]) for (i, t) in enumerate(probes))
    end
end
