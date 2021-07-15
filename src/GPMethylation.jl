module GPMethylation

# Write your package code here.
using CSV: propertynames
using DataFrames, CSV, RData
using GaussianProcesses 
using LinearAlgebra
using ThreadTools
using RCall

export loaddata, loadsamplemeta, loadprobemeta

include("loaddata.jl")

end
