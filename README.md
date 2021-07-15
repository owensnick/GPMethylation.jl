# GPMethylation

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://owensnick.github.io/GPMethylation.jl/dev)


GPMethylation.jl a package that runs Gaussian process regression and model selection on methylation time series.




## Installation
Prior to running the package must be installed using Julia package manager
```julia
    ] add https://github.com/owensnick/GPMethylation.jl
```

## Running
Once installed, example script to run analysis given below:

Make a file `gpanalysis.jl` with the following contents:
```julia
using GPMethylation

samplemetafile = ... # tab separated file describing samples
betafile       = ... # RData file describing probe beta values
outdir         = ... # directory to save results
m              = 50  # number of points at which to calculate GP mean and var between min and max Age

run_gpregression(samplemetafile, betafile, outdir, m)
```

Then at a terminal run
```bash
julia -t N gpanalysis.jl
```
Where `N` is the number of threads.