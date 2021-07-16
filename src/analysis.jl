

"""
    run_gpregression(samplemetafile, betafile, m=50; verbose=true)

Run analysis for model selection between constant, linear or non-linear (Mat52) kernels

Input:
- `samplefmetafile` - File with sample meta
- `betafile` - `RData` for with probe beta.
- `outdir` - Directory to store results, with following contents
  - `gpstats.tsv` - Tab separated file with marginal likelihoods and optimised kernel parameters for all kernels
  - `gpmeanvar.rdata` - RData file containing GP mean and variance for each probe
- `m = 50` - Number of points to calculate GP mean and variance at between min and max Age.
"""
function run_gpregression(samplemetafile, betafile, outdir, m=50; verbose=true)
    nt = Threads.nthreads()
    
    verbose && println("[GPM]\tRunning GP regression on methylation time series with ", nt, " threads.")
    meta, probes, beta = loadall(samplemetafile, betafile, verbose)
    ap = Float64.(meta.Age)
    st = range(minimum(ap), maximum(ap), length=m)
    starttime = time()
    dfgp, μ, v = gpreg_all_threads(ap, beta, st, gpmodels, nt = 16)
    verbose && println("[GPM]\tComplete in ", time() - starttime, " seconds")


    gps = gpstats(probes, dfgp, beta, μ)
    mkpath(outdir)

    f = CSV.write(joinpath(outdir, "gpstats.tsv"), gps, delim='\t')
    verbose && println("[GPM]\tWritten ", f)
    saverdata(probes, st, μ, v, outdir; verbose=verbose)
end


"""
    saverdata(probes, samples, μ, v, outdir; verbose=true)

Save gp mean and var in RData format.
"""
function saverdata(probes, samples, μ, v, outdir; verbose=true)

    file = joinpath(outdir, "gpmeanvar.rdata")
    
    gp_mean = μ'
    gp_var  = v'

    @rput gp_mean gp_var probes samples

    R"""
        colnames(gp_mean) <- samples
        rownames(gp_mean) <- probes
        colnames(gp_var) <- samples
        rownames(gp_var) <- probes

        save(gp_mean, gp_var, file=$file)
    """

    verbose && println("[GPM]\tWritten ", file)
    file

end

"""
    loadresults_rdata(file)

Loads gp_mean and gp_var results from `RData` file back into Julia.
"""
function loadresults_rdata(file)

    R"""
        load($file)
        st <- colnames(gp_mean)
        probes <- rownames(gp_mean)
    """
    @rget gp_mean gp_var st probes

    (μ=Matrix(gp_mean'), v=Matrix(gp_var'), st=parse.(Float64, st), probes=probes)
end


"""
    loadresults(outdir)

Loads gpstats and GP mean and var
"""
function loadresults(outdir)
    gps = CSV.read(joinpath(outdir, "gpstats.tsv"), DataFrame)
    gpres = loadresults_rdata(joinpath(outdir, "gpmeanvar.rdata"))

    gps, gpres
end


"""
    gpstats(probes, dfgp, beta, gpμ)

Writes tab-separated file with gp stats.
"""
function gpstats(probes, dfgp, beta, gpμ)

    @assert length(probes) == size(dfgp, 1)

    probe_min = vec(minimum(beta, 1))
    probe_max = vec(maximum(beta, 1))
    probe_env  = probe_max .- probe_min

    gp_min = vec(minimum(gpμ, 1))
    gp_max = vec(maximum(gpμ, 1))
    gp_env  = gp_max .- gp_min


    evdf = DataFrame(probe_min=probe_min, probe_max=probe_max, probe_env=probe_env,
                     gp_min=gp_min, gp_max=gp_max, gp_env=gp_env)

    [DataFrame(Probe=probes) dfgp evdf]
end


