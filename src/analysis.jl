

"""
    run_gpregression(samplemetafile, betafile, m=50; verbose=true, usepriors=true)

Run analysis for model selection between constant, linear or non-linear (Mat52) kernels

Input:
- `samplefmetafile` - File with sample meta
- `betafile` - `RData` for with probe beta.
- `outdir` - Directory to store results, with following contents
  - `gpstats.tsv` - Tab separated file with marginal likelihoods and optimised kernel parameters for all kernels
  - `gpmeanvar.rdata` - RData file containing GP mean and variance for each probe
- `m = 50` - Number of points to calculate GP mean and variance at between min and max Age.
- usepriors=true - Use default priors on hyperparameters of kernels
"""
function run_gpregression(samplemetafile, betafile, outdir, m=50; verbose=true, usepriors=true)
    nt = Threads.nthreads()
    
    verbose && println("[GPM]\tRunning GP regression on methylation time series with ", nt, " threads.")
    meta, probes, beta = loadall(samplemetafile, betafile, verbose)
    ap = Float64.(meta.Age)
    st = range(minimum(ap), maximum(ap), length=m)
    gpf = ifelse(usepriors, gpmodels_prior, gpmodels)

    if usepriors
        verbose && println("[GPM]\tUsing priors on kernel parameters.")
    end
    starttime = time()
    dfgp, μ, v = gpreg_all_threads(ap, beta, st, gpf, nt = 16)
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

Writes tab-separated file with gp stats. Can accept Union{T, Missing} beta, however, will fail if all values of a column are missing.
"""
function gpstats(probes, dfgp, beta, gpμ)

    @assert length(probes) == size(dfgp, 1)


    probe_min = minimum.(skipmissing.(eachcol(beta)))
    probe_max = maximum.(skipmissing.(eachcol(beta)))
    probe_env  = probe_max .- probe_min

    gp_min = minimum.(skipmissing.(eachcol(gpμ)))
    gp_max = maximum.(skipmissing.(eachcol(gpμ)))
    gp_env  = gp_max .- gp_min


    evdf = DataFrame(probe_min=probe_min, probe_max=probe_max, probe_env=probe_env,
                     gp_min=gp_min, gp_max=gp_max, gp_env=gp_env)

    df = [DataFrame(Probe=probes) dfgp evdf]
     rename!(df, 
        :const_σn    => :const_sign,
        :const_σf    => :const_sigf,
        :linear_mll  => :linear_mll,
        :linear_σn   => :linear_sign,
        :linear_σf   => :linear_sigf,
        :linear_ℓ    => :linear_ell,
        :mat52_mll   => :mat52_mll,
        :mat52_σn    => :mat52_sign,
        :mat52_σf    => :mat52_sigf,
        :mat52_ℓ     => :mat52_ell)    

    df[!, :LR_mat52_const]  = df.mat52_mll .- df.const_mll
    df[!, :LR_mat52_linear] = df.mat52_mll .- df.linear_mll
    df[!, :LR_linear_const] = df.linear_mll .- df.const_mll

    bic_term = log.(vec(sum(.!ismissing.(beta), dims=1)))./2
    df[!, :BIC_mat52_const]  = df.LR_mat52_const - bic_term
    #df[!, :BIC_mat52_linear] = df.LR_mat52_linear ## mat52 and linear models have identical number of params
    df[!, :BIC_linear_const] = df.LR_linear_const - bic_term

    df[!, :PreferedModelLR]   = pref_model.(df.mat52_mll, df.linear_mll, df.const_mll)
    
     df
end

function pref_model(mat52_mll, linear_mll, const_mll, labels=("Mat52", "Linear", "Const"))
    i = argmax((mat52_mll, linear_mll, const_mll))
    labels[i]
end


