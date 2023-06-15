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


noise_labels = ["damping", "no_noise"]
noiseless_return = "Noiseless"
gate_labels = ["PauliX", "HGate", "HCZGate"]
gateless_return = "Gateless"
num_qubits = ["Single","Two"]
qubitless_return = "Qubitless"

paths = (data = data_path,figs = fig_path)

p = (noise_labels = noise_labels,
    noiseless_return = noiseless_return,
    gate_labels = gate_labels,
    gateless_return = gateless_return,
    num_qubits = num_qubits,
    qubitless_return = qubitless_return)


loaded_data_saved_plots = get_save_probability_plots(paths,p)

loaded_data_saved_plots


function get_save_probability_plots(paths,p)
    data_path,fig_path = paths
    noise_labels, noiseless_return, gate_labels, gateless_return, num_qubits, qubitless_return = p
    df = @chain data_path begin
        find(_,0,4,".csv")    
        filter(contains("probabilities"),_) 
        DataFrame(path = _) 
        @rtransform :basename = basename(:path) 
        @rtransform :num_qubit_regex = rjoin_or(num_qubits) 
        @rtransform :qubit_label = get_label(:num_qubit_regex,:basename,qubitless_return) 
        @rtransform :noise_regex = rjoin_or(noise_labels)
        @rtransform :noise_label = get_label(:noise_regex,:basename,noiseless_return) 
        @rtransform :gate_regex = rjoin_or(gate_labels) 
        @rtransform :gate_label = get_label(:gate_regex,:basename,gateless_return) 
        @rtransform :noise_prob = get_noise_prob(:basename) 
        @rtransform :prob_df = CSV.read(:path,DataFrame,header=true)
        @rtransform :x_axis = 1:size(:prob_df,1)
        @rtransform :num_qubits = Int(size(:prob_df,1)/2)
        @rtransform :enum_comp_basis = list_qubit_axis_arrangement(:num_qubits)
        @rtransform :probs = :prob_df.probabilities 
        @rtransform :plots = plot_probability_barplot(:x_axis,:probs,:enum_comp_basis,:noise_prob,:gate_label,:noise_prob)    
        @rtransform :path1 = join([:noise_label,:qubit_label,"Qubits",:gate_label])
        @rtransform :path2 = join(["_NoiseModel-",:noise_prob])
        @rtransform :plot_file_name = join([:path1,:path2,".png"])
        @select $(Not([:path1,:path2]))
        @rtransform :plot_str = joinpath(fig_path,:plot_file_name)
        @rtransform :plot_file_res = save_plot_verify(:plot_str,:plots)
    end
    return df
end

    

using SparseArrays

# Define the size of the matrix
n = 4  # Number of qubits

# Define the indices and values for the non-zero elements of the matrix
rows = [1, 2, 3, 4]
cols = [1, 2, 3, 4]
vals = [1, 1, 1, -1]

# Create the sparse matrix
cz_matrix = sparse(rows, cols, vals, n, n)
ket0 = [1,0]
H = [1 1;1 -1]
pstate = H*ket0
pstate'
istate = kron(pstate,pstate)

pstate'*cz_matrix*istate