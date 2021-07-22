module GPMethylation

using Statistics
using DataFrames, CSV, RData
using GaussianProcesses
using Distributions 
using LinearAlgebra
using ThreadTools
using RCall
using ProgressMeter

export  loaddata, loadsamplemeta, loadprobemeta, loadall, loadmeta,
        gpreg, gpreg_const, gpreg_linear, gpreg_matern52, gpreg_all_threads, gpmodels,
        run_gpregression, saverdata, gpstats, loadresults
       

include("loaddata.jl")
include("gp.jl")
include("analysis.jl")

end
