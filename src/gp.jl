

"""
    gpreg_all_threads(meta, beta, gp_reg, ; nt = 16, bt = 1, n = size(R, 1), gpargs...)

Wrapper to run a Gaussian Process regression function `gp_reg ∈ {gp_reg_matern52, gp_reg_const, gp_reg_linear}` with multiple threads.
    
Keyword argument `nt = 16` specifies number of threads. 

Allows balance of threads to run regression on independent probes in parallel vs threads used by BLAS library for parallelising linear algebra used in regression of a single probe.

Currently, recommend setting BLAS threads `bt = 1` and parallelising over multiple probes.

"""
function gpreg_all_threads(meta, beta, gp_reg, ; nt = 16, bt = 1, n = size(R, 1), gpargs...)
        BLAS.set_num_threads(bt)
        ap = Float64.(meta.Age)
        st = range(minimum(meta.Age), maximum(meta.Age), length=50)
        GPT = tmap(i -> gp_reg(ap, beta[i, :], st, gpargs...), nt, 1:n)
        BLAS.set_num_threads(btc)
        GPT
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
Named tuple (gp=gp, μ=μ, v=v, st=st)
"""
function gpreg(t, y, st, kernel, logσn=log(std(y)))
    gp = GP(t, y, MeanZero(), kernel, logσn)
    optimize!(gp)
    μ, v = predict_y(gp, st)
    (gp=gp, μ=μ, v=v, st=st)
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
function gpreg_const(t, y, st, kernparams=(1.0, 1.0), logσn=log(std(y)))
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

