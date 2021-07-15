

"""
    gpreg_all_threads(t, beta, st, gpreg, ; nt = 16, bt = 1, n = size(beta, 2))

Wrapper to run a Gaussian Process regression function `gp_reg ∈ {gp_reg_matern52, gp_reg_const, gp_reg_linear}` with multiple threads.
    
Keyword argument `nt = 16` specifies number of threads. 

Allows balance of threads to run regression on independent probes in parallel vs threads used by BLAS library for parallelising linear algebra used in regression of a single probe.

Currently, recommend setting BLAS threads `bt = 1` and parallelising over multiple probes.

"""
function gpreg_all_threads(t, beta, st, gpreg, ; nt = 16, bt = 1, n = size(beta, 2))
        BLAS.set_num_threads(bt)
        GPT = tmap(i -> gpreg(t, beta[:, i], st), nt, 1:n)
        BLAS.set_num_threads(btc)
        GPT
end


"""
    gpmodels(t, y, st)

Runs GP regression for constant, linear and Matern52 kernels.
"""
function gpmodels(t, y, st)

    ### for now hardcode const, linear and mat52 regression
    ### improvement would be to specify vector of kernels and loop over

    gpc = gpreg_const(t, y, st)
    gpl = gpreg_linear(t, y, st)
    gpm = gpreg_matern52(t, y, st)

    params = (const_mll    = gpc.gp.mll,
     const_σn2    = exp(2*gpc.gp.logNoise.value),
     const_σf2    = gpc.gp.kernel.σ2,
     linear_mll   = gpl.gp.mll,
     linear_σn2   = exp(2*gpl.gp.logNoise.value),
     linear_σf2   = gpl.gp.kernel.σ2,
     linear_ℓ     = gpl.gp.kernel.ℓ,
     mat52_mll    = gpm.gp.mll,
     mat52_σn2    = exp(2*gpm.gp.logNoise.value),
     mat52_σf2    = gpm.gp.kernel.σ2,
     mat52_ℓ      = gpm.gp.kernel.ℓ)

     gpmeanvar = (mat52_μ   = gpm.μ, mat52_v  = gpm.v)

     params, gpmeanvar

end


"""
    gpreg(t, y, st, kernel, kernparams=(1.0, 1.0), logNoise=log(std(y)))

Generic GP stationary zero mean regression with Gaussian likelihood. Params:
  - `t` vector of time points (in this case ages)
  - `y` vector of signal values (in this case prob methylation)
  - `st` vector of time points to evaluate GP `mean` and `var`
  - `kernel` GP `kernel` with inital params e.g. `Mat52Iso`
  - `logσn` log(σn)

Returns:
Named tuple (gp=gp, μ=μ, v=v)
"""
function gpreg(t, y, st, kernel, logσn=log(std(y)))
    gp = GP(t, y, MeanZero(), kernel, logσn)
    optimize!(gp)
    μ, v = predict_y(gp, st)
    (gp=gp, μ=μ, v=v)
end

"""
    gpreg_matern52(t, y, st, kernparams=(1.0, 1.0) , log(std(y))))

GP regression with `Mat52Iso(ℓ, σf)` kernel.

`kernparms = (logℓ, log(σf))`

`ℓ`  - initial lengthscale

`σf` - initial function variance
"""
function gpreg_matern52(t, y, st, kernparams=(1.0, 1.0), logσn=log(std(y)))
    logℓ, logσf = kernparams
    gpreg(t, y, st, Mat52Iso(logℓ, logσf), logσn)
end

"""
    gpreg_const(t, y, st, kernparams=(1.0) , log(std(y))))

GP regression with `Const(σf)`.

`kernparms = (log(σf))`

`σf` - initial function variance
"""
function gpreg_const(t, y, st, kernparams=(1.0), logσn=log(std(y)))
    logσf = kernparams
    gpreg(t, y, st, Const(logσf), logσn)
end  

"""
    gpreg_linear(t, y, st, kernparams=(1.0) , log(std(y))))

GP regression with `Const(σf) + LinIso(ℓ)` kernel.

`kernparms = (logℓ, log(σf))`

`ℓ`  - initial lengthscale

`σf` - initial function variance
"""
function gpreg_linear(t, y, st, kernparams=(1.0, 1.0), logσn=log(std(y)))
    logℓ, logσf = kernparams
    gpreg(t, y, st, Const(logσf) + LinIso(logℓ), logσn)
end  

