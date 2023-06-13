using Pkg
Pkg.activate(".")
using Find
using LinearAlgebra
using CSV
using DelimitedFiles
using DataFrames
using DataFramesMeta
using Chain
using CairoMakie
using BenchmarkTools
using Formatting
using Distributions
using Random
attributes = Theme(fontsize = 25,resolution=(800,800),Axis = (; aspect = 1),Legend = (;framevisible=false))
set_theme!(attributes)


include("functions.jl")


data_path = string(homedir(),"/Projects/exploringNoiseInQuEST/data")
fig_path = string(homedir(),"/Projects/exploringNoiseInQuEST/figs")



df = @chain data_path begin
    find(_,0,4,".csv")    
    filter(contains(r"damping|no_noise"),_) 
    filter(contains("Single"),_) 
    filter(contains("probabilities"),_) 
    DataFrame(path = _)
    @rtransform :basename = basename(:path)
    @rtransform :noise_prob = get_noise_prob(:basename)
    @rtransform :prob_df = CSV.read(:path,DataFrame,header=true)
    @rtransform :x_axis = 1:size(:prob_df,1)
    @rtransform :num_qubits = Int(size(:prob_df,1)/2)
    @rtransform :enum_comp_basis = list_qubit_axis_arrangement(:num_qubits)
    @rtransform :probs = :prob_df.probabilities
    @rtransform :plots = plot_probability_barplot(:x_axis,:probs,:enum_comp_basis,:noise_prob,:basename)    
end
