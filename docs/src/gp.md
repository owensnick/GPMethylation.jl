# Gaussian Process Regression

```@docs
gpreg_all_threads(t, beta, st, gpreg=gpmodels ; nt = 16, bt = 1, n = size(beta, 2))
gpmodels(t, y, st)
gpmodels_prior(td, yd, st)
gpreg(t, y, st, kernel, logσn=log(std(y)))
gpreg_matern52(t, y, st, kernparams=(1.0, 1.0), logσn=log(std(y)))
gpreg_matern52_prior(t, y, st, kernparams=(2.0, log(mean(y))), logσn=log(std(y)), kernpriors=[Normal(4.7, 1.2), Normal(-1, 6)])
gpreg_const(t, y, st, kernparams=(1.0), logσn=log(std(y)))
gpreg_linear(t, y, st, kernparams=(1.0, 1.0), logσn=log(std(y)))
gpreg_linear_prior(t, y, st, kernparams=(1.0, log(mean(y))), logσn=log(std(y)), kernpriors=[Normal(5.7, 2.2)])
```