# The netcdf reader's open-file cache. Like the Provider's memo it may change how
# often a file is OPENED and nothing else, so reads through a kept handle are
# compared with `isequal` on every array against reads that open the file afresh.
# Reuses the fixtures and helpers of test_provider_cache.jl (included first).

nc_opens() = EarthSciIO._NC_OPENS[]
nc_live() = EarthSciIO._NC_LIVE[]
held() = EarthSciIO._nc_held()
kept(path) = EarthSciIO._NC_HANDLES[path]
handles_off(f) = withenv(f, "EARTHSCIIO_NC_HANDLES" => "0")
nc_isopen(h) = isopen(h.ds)

@testset "netcdf reader — open-file cache" begin
    EarthSciIO._nc_close_all!()                  # whatever earlier testsets left open
    @test isempty(held()) && nc_live() == 0
    c, root = stage_days()
    paths = [blob_path(root, u) for u in DAY_URLS]
    R = FORMAT_REGISTRY["netcdf"]
    SEL = Dict("axes" => Any["all", Dict("indices" => [2, 0]), Dict("slice" => [1, 3])])

    @testset "one open per file, however many reads" begin
        n = nc_opens()
        for _ in 1:5, p in paths
            read_native(R, p; variables = ["t2m"])
        end
        @test nc_opens() == n + 3
        @test held() == paths
        @test nc_live() == 3
        @test dim_length(R, paths[1], "time") == 4
        @test dim_length(R, paths[1], "nope") === nothing
        @test nc_opens() == n + 3                # dim_length rides the same handle
    end

    @testset "a kept handle reads what a fresh open reads" begin
        rec = Dict{String,Any}("dim" => "time", "indices" => [3, 0])
        cases = [NamedTuple(), (variables = ["t2m"],), (variables = String[],),
                 (select = SEL,), (records = rec,),
                 (records = Dict{String,Any}("dim" => "time", "indices" => len -> [len - 1]),),
                 (variables = ["t2m"], select = SEL, records = rec)]
        # The references first: a capacity of 0 closes whatever is kept.
        want = handles_off() do
            Dict((i, j) => read_native(R, paths[i]; cases[j]...)
                 for i in eachindex(paths), j in eachindex(cases))
        end
        n = nc_opens()
        for _ in 1:3, i in eachindex(paths), j in eachindex(cases)
            @test identical_native(read_native(R, paths[i]; cases[j]...), want[(i, j)])
        end
        @test nc_opens() == n + 3                # one open per file; the rest were kept
    end

    @testset "Provider samples: every combination of the two caches == neither" begin
        probes = [0.0, 0.5, 3.0, 3.0, 4.0, 7.5, 11.0, 30.0, 11.0, 8.0, 7.0, 3.5, 0.0]
        combos = [("1", "32"), ("0", "32"), ("1", "0"), ("1", "1")]
        for rps in (nothing, 1, 2), sel in (nothing, SEL), push in (true, false)
            kw = push ? NamedTuple() : (reader_kwargs = (records = nothing,),)
            p = discrete_provider(c, day_url, DAY_TICKS; format = "netcdf",
                                  time_dim = "time", records_per_sample = rps, kw...)
            for t in probes
                want = withenv("EARTHSCIIO_PROVIDER_CACHE" => "0",
                               "EARTHSCIIO_NC_HANDLES" => "0") do
                    materialize(p, t; select = sel)
                end
                for (memo, nh) in combos
                    got = withenv("EARTHSCIIO_PROVIDER_CACHE" => memo,
                                  "EARTHSCIIO_NC_HANDLES" => nh) do
                        materialize(p, t; select = sel)
                    end
                    @test identical_native(got, want)
                end
            end
        end
    end

    @testset "a touched or rewritten file is reopened, and the old handle closed" begin
        _, root2 = stage_days()
        p = blob_path(root2, DAY_URLS[1])
        cell() = read_native(R, p; variables = ["t2m"])["t2m"].data[:, 1, 1]
        @test cell() ≈ [101.11, 102.11, 103.11, 104.11]
        n = nc_opens()
        h = kept(p)
        @test cell() ≈ [101.11, 102.11, 103.11, 104.11]
        @test nc_opens() == n && kept(p) === h   # warm

        backdate!(p, 30)                         # touched: same bytes, new mtime
        @test cell() ≈ [101.11, 102.11, 103.11, 104.11]
        @test nc_opens() == n + 1 && !nc_isopen(h)
        h = kept(p)

        tmp = write_hourly_nc(tempname() * ".nc"; nrec = 4, base = 500.0, t0 = 0.0)
        @test filesize(tmp) == filesize(p)
        mv(tmp, p; force = true)                 # replaced by a rename, same size
        backdate!(p, 45)
        @test cell() ≈ [501.11, 502.11, 503.11, 504.11]
        @test nc_opens() == n + 2 && !nc_isopen(h)
        h = kept(p)

        tmp = write_hourly_nc(tempname() * ".nc"; nrec = 4, base = 700.0, t0 = 0.0)
        write(p, read(tmp))                      # rewritten in place: same inode, size
        backdate!(p, 90)
        @test cell() ≈ [701.11, 702.11, 703.11, 704.11]
        @test nc_opens() == n + 3 && !nc_isopen(h)

        make_racy!(p)                            # inside the margin: never kept
        cell(); cell()
        @test nc_opens() == n + 5
        @test !(p in held())
    end

    @testset "touching or rewriting a blob invalidates both caches" begin
        c3, root3 = stage_days()
        path = blob_path(root3, DAY_URLS[1])
        p = discrete_provider(c3, day_url, DAY_TICKS; format = "netcdf-counting",
                              time_dim = "time", records_per_sample = 2)
        cell() = materialize(p, 1.0)["t2m"].data[:, 1, 1]
        @test cell() ≈ [102.11, 103.11]
        d, o = decodes(), nc_opens()
        @test cell() ≈ [102.11, 103.11]
        @test (decodes(), nc_opens()) == (d, o)            # both warm
        backdate!(path, 30)
        @test cell() ≈ [102.11, 103.11]
        @test (decodes(), nc_opens()) == (d + 1, o + 1)
        mv(write_hourly_nc(tempname() * ".nc"; nrec = 4, base = 500.0, t0 = 0.0), path;
           force = true)
        backdate!(path, 45)
        @test cell() ≈ [502.11, 503.11]
        @test (decodes(), nc_opens()) == (d + 2, o + 2)
        write(path, read(write_hourly_nc(tempname() * ".nc"; nrec = 4, base = 700.0, t0 = 0.0)))
        backdate!(path, 90)
        @test cell() ≈ [702.11, 703.11]
        @test (decodes(), nc_opens()) == (d + 3, o + 3)
    end

    @testset "LRU eviction closes the least recently used handle" begin
        ref2 = handles_off(() -> read_native(R, paths[2]))
        withenv("EARTHSCIIO_NC_HANDLES" => "2") do
            read_native(R, paths[1]); read_native(R, paths[2])
            h1, h2 = kept(paths[1]), kept(paths[2])
            read_native(R, paths[1])             # paths[1] is now the most recent
            read_native(R, paths[3])             # ...so paths[2] goes
            @test held() == [paths[1], paths[3]]
            @test !nc_isopen(h2) && nc_isopen(h1)
            @test nc_live() == 2
            n = nc_opens()
            @test identical_native(read_native(R, paths[2]), ref2)
            @test nc_opens() == n + 1            # reopened, and paths[1] went
            @test held() == [paths[3], paths[2]] && !nc_isopen(h1)
        end
        withenv("EARTHSCIIO_NC_HANDLES" => "1") do   # a smaller capacity evicts on the next read
            read_native(R, paths[2])
            @test held() == [paths[2]] && nc_live() == 1
        end
    end

    @testset "capacity 0 opens per read, and closes what was kept" begin
        read_native(R, paths[1])
        h = kept(paths[1])
        handles_off() do
            n = nc_opens()
            a = read_native(R, paths[1])
            b = read_native(R, paths[1])
            @test nc_opens() == n + 2
            @test isempty(held()) && nc_live() == 0
            @test !nc_isopen(h)
            @test identical_native(a, b)
        end
        for bad in ("-1", "lots", "3.5")
            @test_throws ArgumentError withenv(() -> read_native(R, paths[1]),
                                               "EARTHSCIIO_NC_HANDLES" => bad)
        end
    end

    @testset "only cache blobs are kept, so any other file can still be rewritten" begin
        p = backdate!(write_hourly_nc(tempname() * ".nc"; nrec = 4, base = 100.0, t0 = 0.0))
        n = nc_opens()
        read_native(R, p); read_native(R, p)
        @test nc_opens() == n + 2 && !(p in held())
        # HDF5 refuses to reopen for writing a netCDF-4 file this process holds
        # open; nothing is held, so the rewrite goes through.
        write_hourly_nc(p; nrec = 4, base = 900.0, t0 = 0.0)
        @test read_native(R, p)["t2m"].data[1, 1, 1] ≈ 901.11
        @test EarthSciIO._nc_keepable(paths[1])
        @test !EarthSciIO._nc_keepable(p)
        @test !EarthSciIO._nc_keepable(joinpath(dirname(paths[1]), "notakey.nc"))
        k = cache_key("https://elsewhere")       # a key under the wrong shard
        @test !EarthSciIO._nc_keepable(joinpath(root, "v1", "blobs",
                                                k[1:2] == "00" ? "01" : "00", k * ".nc"))
    end

    @testset "a failed read does not keep its handle" begin
        read_native(R, paths[1])
        @test paths[1] in held()
        @test_throws ArgumentError read_native(R, paths[1]; variables = ["nope"])
        @test !(paths[1] in held())
        n = nc_opens()
        read_native(R, paths[1])
        @test nc_opens() == n + 1 && paths[1] in held()
    end

    @testset "concurrent reads of one file, with and without eviction under them" begin
        want = Dict(p => handles_off(() -> read_native(R, p; variables = ["t2m"]))
                    for p in paths)
        for cap in ("32", "1")
            withenv("EARTHSCIIO_NC_HANDLES" => cap) do
                jobs = repeat(paths, 16)
                got = Vector{Any}(undef, length(jobs))
                @sync for (i, p) in enumerate(jobs)
                    Threads.@spawn got[i] = read_native(R, p; variables = ["t2m"])
                end
                @test all(identical_native(got[i], want[p]) for (i, p) in enumerate(jobs))
                @test length(held()) <= parse(Int, cap)
                @test nc_live() == length(held())    # every handle that left was closed
            end
        end
    end

    @testset "handles are closed at exit" begin
        @test EarthSciIO._nc_close_all! in Base.atexit_hooks
        EarthSciIO._nc_close_all!()
        @test isempty(held()) && nc_live() == 0
    end
end
