

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
    GPT = gpreg_all_threads(ap, beta, st, gpmodels, nt = 16)
    verbose && println("[GPM]\tComplete in ", time() - starttime, " seconds")


    gps = gpstats(probes, GPT, beta)
    mkpath(outdir)

    f = CSV.write(joinpath(outdir, "gpstats.tsv"), gps, delim='\t')
    verbose && println("[GPM]\tWritten ", f)
    saverdata(probes, st, GPT, outdir; verbose=verbose)
end


"""
    saverdata(probes, samples, GPT, outdir; verbose=true)

Save gp mean and var in RData format.
"""
function saverdata(probes, samples, GPT, outdir; verbose=true)

    file = joinpath(outdir, "gpmeanvar.rdata")
    

    gp_mean = mapreduce(g -> g[2].mat52_μ, hcat, GPT)';
    gp_var  = mapreduce(g -> g[2].mat52_v, hcat, GPT)';

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

    (μ=gp_mean, v=gp_var, st=st, probes=probs)
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
    gpstats(probes, GPT, beta)  

Writes tab-separated file with gp stats.
"""
function gpstats(probes, GPT, beta)

    @assert length(probes) == length(GPT)

    gpdf = DataFrame(
     const_mll    = Float64[],
     const_σn2    = Float64[],
     const_σf2    = Float64[],
     linear_mll   = Float64[],
     linear_σn2   = Float64[],
     linear_σf2   = Float64[],
     linear_ℓ     = Float64[],
     mat52_mll    = Float64[],
     mat52_σn2    = Float64[],
     mat52_σf2    = Float64[],
     mat52_ℓ      = Float64[])

     evdf = DataFrame(
     probe_min     = Float64[],
     probe_max     = Float64[],
     probe_env     = Float64[],
     gp_min        = Float64[],
     gp_max        = Float64[],
     gp_env        = Float64[])


    for i = 1:length(probes)

        gpparams = GPT[i][1]
        gp_meanvar = GPT[i][2]

        probe_min = minimum(beta[:, i])
        probe_max = maximum(beta[:, i])
        probe_ev  = probe_max - probe_min

        gp_min = minimum(gp_meanvar.mat52_μ)
        gp_max = maximum(gp_meanvar.mat52_μ)
        gp_ev  = gp_max - gp_min

        
        push!(gpdf, gpparams)
        push!(evdf, (probe_min, probe_max, probe_ev, gp_min, gp_max, gp_ev))
    end

    [DataFrame(Probe=probes) gpdf evdf]
end


