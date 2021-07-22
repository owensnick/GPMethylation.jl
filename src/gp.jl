

"""
gpreg_all_threads(t, beta, st, gpreg=gpmodels ; nt = 16, bt = 1, n = size(beta, 2))

Wrapper to run a Gaussian Process regression function `gp_reg ∈ {gp_reg_matern52, gp_reg_const, gp_reg_linear}` with multiple threads.
    
Allows balance of threads to run regression on independent probes in parallel vs threads used by BLAS library for parallelising linear algebra used in regression of a single probe.

Currently, recommend setting BLAS threads `bt = 1` and parallelising over multiple probes.

"""
function gpreg_all_threads(t, beta, st, gpreg=gpmodels ; nt = 16, bt = 1, n = size(beta, 2))
        btc = ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
        BLAS.set_num_threads(bt)

    dft = DataFrame(const_mll    = Vector{Float64}(undef, n),
                    const_σn2    = Vector{Float64}(undef, n),
                    const_σf2    = Vector{Float64}(undef, n),
                    linear_mll   = Vector{Float64}(undef, n),
                    linear_σn2   = Vector{Float64}(undef, n),
                    linear_σf2   = Vector{Float64}(undef, n),
                    linear_ℓ     = Vector{Float64}(undef, n),
                    mat52_mll    = Vector{Float64}(undef, n),
                    mat52_σn2    = Vector{Float64}(undef, n),
                    mat52_σf2    = Vector{Float64}(undef, n),
                    mat52_ℓ      = Vector{Float64}(undef, n))
    μ = Matrix{Float64}(undef, length(st), n)
    v = Matrix{Float64}(undef, length(st), n)

    p = Progress(n, barlen=75)
    Threads.@threads for i = 1:n
        params, gpu, gpv = gpreg(t, beta[:, i], st)
        dft[i, :] = params
        μ[:, i] .= gpu
        v[:, i] .= gpv
        next!(p)
    end

    BLAS.set_num_threads(btc)
    dft, μ, v
end

filtmissing(t, y) = t, y
function filtmissing(t, y::Vector{Union{T, Missing}}) where {T}
    ind = .!ismissing.(y)
    t[ind], Vector{T}(y[ind])
end



"""
    gpmodels(t, y, st)

Runs GP regression for constant, linear and Matern52 kernels.
"""
function gpmodels(td, yd, st)

    ### for now hardcode const, linear and mat52 regression
    ### improvement would be to specify vector of kernels and loop over

    ### remove missings if necessary
    t, y = filtmissing(td, yd)

    gpc = gpreg_const(t, y, st)
    gpl = gpreg_linear(t, y, st)
    gpm = gpreg_matern52(t, y, st)

    params = (const_mll    = gpc.gp.mll,
     const_σn    = exp(gpc.gp.logNoise.value),
     const_σf    = sqrt(gpc.gp.kernel.σ2),
     linear_mll  = gpl.gp.mll,
     linear_σn   = exp(gpl.gp.logNoise.value),
     linear_σf   = sqrt(gpl.gp.kernel.kleft.σ2),
     linear_ℓ    = sqrt(gpl.gp.kernel.kright.ℓ2),
     mat52_mll   = gpm.gp.mll,
     mat52_σn    = exp(gpm.gp.logNoise.value),
     mat52_σf    = sqrt(gpm.gp.kernel.σ2),
     mat52_ℓ     = gpm.gp.kernel.ℓ)

     params, gpm.μ,  gpm.v
end

"""
    gpmodels_prior(t, y, st)

Runs GP regression for constant, linear and Matern52 kernels with priors on hyperparameters for constant and linear kernels
"""
function gpmodels_prior(td, yd, st)

    ### for now hardcode const, linear and mat52 regression
    ### improvement would be to specify vector of kernels and loop over

    ### remove missings if necessary
    t, y = filtmissing(td, yd)

    gpc = gpreg_const(t, y, st)
    gpl = gpreg_linear(t, y, st)
    gpm = gpreg_matern52_prior(t, y, st)

    params = (const_mll    = gpc.gp.mll,
     const_σn    = exp(gpc.gp.logNoise.value),
     const_σf    = sqrt(gpc.gp.kernel.σ2),
     linear_mll  = gpl.gp.mll,
     linear_σn   = exp(gpl.gp.logNoise.value),
     linear_σf   = sqrt(gpl.gp.kernel.kleft.σ2),
     linear_ℓ    = sqrt(gpl.gp.kernel.kright.ℓ2),
     mat52_mll   = gpm.gp.mll,
     mat52_σn    = exp(gpm.gp.logNoise.value),
     mat52_σf    = sqrt(gpm.gp.kernel.σ2),
     mat52_ℓ     = gpm.gp.kernel.ℓ)

     params, gpm.μ,  gpm.v
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
    gpreg_matern52(t, y, st, kernparams=(1.0, log(mean(y))) , log(std(y))))

GP regression with `Mat52Iso(ℓ, σf)` kernel.

`kernparms = (logℓ, log(σf))`

`ℓ`  - initial lengthscale

`σf` - initial function variance
"""
function gpreg_matern52_prior(t, y, st, kernparams=(2.0, log(mean(y))), logσn=log(std(y)), kernpriors=[Normal(4.7, 1.2), Normal(-1, 6)])
    logℓ, logσf = kernparams
    kern = Mat52Iso(logℓ, logσf)
    set_priors!(kern, kernpriors)
    gpreg(t, y, st, kern, logσn)
end

"""
    gpreg_matern52(t, y, st, kernparams=(1.0, 1.0) , log(std(y))))

GP regression with `Mat52Iso(ℓ, σf)` kernel.

`kernparms = (logℓ, log(σf))`

`ℓ`  - initial lengthscale

`σf` - initial function variance
"""
function gpreg_matern52(t, y, st, kernparams=(1.0, log(mean(y))), logσn=log(std(y)))
    logℓ, logσf = kernparams
    gpreg(t, y, st, Mat52Iso(logℓ, logσf), logσn)
end

"""
    gpreg_const(t, y, st, kernparams=(1.0) , log(std(y))))

GP regression with `Const(σf)`.

`kernparms = (log(σf))`

`σf` - initial function variance
"""
function gpreg_const(t, y, st, kernparams=(log(mean(y))), logσn=log(std(y)))
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
function gpreg_linear(t, y, st, kernparams=(1.0, log(mean(y))), logσn=log(std(y)))
    logℓ, logσf = kernparams
    gpreg(t, y, st, Const(logσf) + LinIso(logℓ), logσn)
end  

"""
    gpreg_linear_prior(t, y, st, kernparams=(1.0) , log(std(y))))

GP regression with `Const(σf) + LinIso(ℓ)` kernel.

`kernparms = (logℓ, log(σf))`

`ℓ`  - initial lengthscale

`σf` - initial function variance
"""
function gpreg_linear_prior(t, y, st, kernparams=(1.0, log(mean(y))), logσn=log(std(y)), kernpriors=[Normal(5.7, 2.2)])
    logℓ, logσf = kernparams
    kc = Const(logσf)
    kl = LinIso(logℓ)
    set_priors!(kl, kernpriors)
    k = kc + kl
    gpreg(t, y, st, k, logσn)
end  