



function gp_reg_all_threads(meta, R, gp_reg, ; nt = 16, bt = 1, n = size(R, 1), gpargs...)
        BLAS.set_num_threads(bt)
        ap = Float64.(meta.Age)
        st = range(minimum(meta.Age), maximum(meta.Age), length=50)
        GPT = tmap(i -> gp_reg(ap, R[i, :], st, gpargs...), nt, 1:n)
        BLAS.set_num_threads(btc)
        GPT
end


function gp_reg_matern52(t, y, st, initparams=(1.0, 1.0, log(std(y))))
    σf, ℓ, σn = initparams
    kern = Mat52Iso(ℓ, σf)
    gp = GP(t, y, MeanZero(), kern, σn)
    optimize!(gp)
    μ, v = predict_y(gp, st)
    (gp=gp, μ=μ, v=v, st=st)
end


function gp_reg_const(t, y, st; initparams = (1.0, 1.0))
    σf, σn = initparams
    kern = Const(σf)
    gp = GP(t, y, MeanZero(), kern, σn)
    optimize!(gp)
    μ, v = predict_y(gp, st)
    (gp=gp, μ=μ, v=v, st=st)
end



function gp_reg_linear(t, y, st; initparams = (1.0, 1.0, 1.0))
    σf, ℓ, σn = initparams
    kern = Const(σf) + LinIso(ℓ)
    gp = GP(t, y, MeanZero(), kern, σn)
    optimize!(gp)
    μ, v = predict_y(gp, st)
    (gp=gp, μ=μ, v=v, st=st)
end