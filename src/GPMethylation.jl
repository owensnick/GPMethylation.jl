module GPMethylation

# Write your package code here.
using CSV: propertynames
using DataFrames, CSV, RData
using GaussianProcesses 
using LinearAlgebra
using ThreadTools
using RCall

export loaddata, loadsamplemeta, loadprobemeta, gpreg, gpreg_const, gpreg_linear, gpreg_matern52, gpreg_all_threads

include("loaddata.jl")
include("gp.jl")

end
