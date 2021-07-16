# Gaussian Process Regression

```@docs
gpreg_all_threads(t, beta, st, gpreg=gpmodels ; nt = 16, bt = 1, n = size(beta, 2))
gpmodels(t, y, st)
gpreg(t, y, st, kernel, logσn=log(std(y)))
gpreg_matern52(t, y, st, kernparams=(1.0, 1.0), logσn=log(std(y)))
gpreg_const(t, y, st, kernparams=(1.0), logσn=log(std(y)))
gpreg_linear(t, y, st, kernparams=(1.0, 1.0), logσn=log(std(y)))
```