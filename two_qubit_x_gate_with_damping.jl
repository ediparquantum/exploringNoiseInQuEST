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


# Example usage

single_qubit_axis_vec = list_qubit_axis_arrangement(1)
two_qubit_axis_vec = list_qubit_axis_arrangement(2)




@chain data_path begin
    find(_,0,4,".csv")    
    @aside @chain _ begin
        one_qubit = filter(contains("single"),_) 
    end
    @aside @chain _ begin
        two_qubit = filter(contains("two"),_)
    end
end

"""
    Plot single qubit probabilities as bar plots
"""
prob_path_single = filter(contains("probabilities"),one_qubit)
xgate_single_df = CSV.read(prob_path_single[1],DataFrame,header=true)
xgate_damp_single_df = CSV.read(prob_path_single[2],DataFrame,header=true)

x_axis = 1:size(xgate_single_df,1)
height = xgate_single_df.probabilities
f=Figure()
ax=Axis(f[1,1],xlabel=L"|\psi\rangle",ylabel = L"P(|\psi\rangle)",xticks = (x_axis, single_qubit_axis_vec))
barplot!(ax,x_axis,height,width=0.4,label="Pauli X")
ylims!(0,1)
axislegend("Single Qubit",position=:lt)
f
save(fig_path*"/single_qubit/x-gate-probs.png",f)


x_axis = 1:size(xgate_single_df,1)
height = xgate_damp_single_df.probabilities
f=Figure()
ax=Axis(f[1,1],xlabel=L"|\psi\rangle",ylabel = L"P(|\psi\rangle)",xticks = (x_axis, single_qubit_axis_vec))
barplot!(ax,x_axis,height,width=0.4,label="Pauli X\nDamping, p=0.1")
ylims!(0,1)
axislegend("Single Qubit",position=:lt)
f
save(fig_path*"/single_qubit/x-gate-damp-0.1-probs.png",f)



"""
    Plot two qubit probabilities as bar plots
"""
prob_path_two = filter(contains("probabilities"),two_qubit)
xgate_two_df = CSV.read(prob_path_two[1],DataFrame,header=true)
xgate_damp_two_df = CSV.read(prob_path_two[2],DataFrame,header=true)


x_axis = 1:size(xgate_two_df,1)
height = xgate_two_df.probabilities
f=Figure()
ax=Axis(f[1,1],xlabel=L"|\psi\rangle",ylabel = L"P(|\psi\rangle)",xticks = (x_axis, two_qubit_axis_vec))
barplot!(ax,x_axis,height,width=0.4,label="Pauli X")
ylims!(0,1)
axislegend("Two Qubits",position=:lt)
f
save(fig_path*"/two_qubit/x-gate-probs.png",f)

x_axis = 1:size(xgate_two_df,1)
height = xgate_damp_two_df.probabilities
f=Figure()
ax=Axis(f[1,1],xlabel=L"|\psi\rangle",ylabel = L"P(|\psi\rangle)",xticks = (x_axis, two_qubit_axis_vec))
barplot!(ax,x_axis,height,width=0.4,label="Pauli X\nDamping, p=0.1")
ylims!(0,1)
axislegend("Two Qubits",position=:lt)
f



"""
    Test X gate density matrix probabilities
"""
# Simulated in QuEST
x_gate = "data/amplitudes_single_qubit-X.csv"
x_gate_quest_probs = "data/probabilities_single_qubit-X.csv"
quest_probs = CSV.read(x_gate_quest_probs,DataFrame,header=true)
amp_mat = load_amp_vector_density_mat(x_gate)
prob_vec = get_prob_vec(amp_mat) 



# Simulated in Julia
ρ = Complex.([1,0])*transpose(Complex.([1,0]))
X = [0 1;1 0]
ρ̂ = update_ρ(ρ,X)
prob_vec_sim = get_prob_vec(ρ̂)

@chain quest_probs begin
    @rename :QuEST = $1
    @transform :JuliaQuEST = prob_vec
    @transform :Julia = prob_vec_sim
end

"""
    Test X gate density matrix probabilities with damping
"""
# Simulated in QuEST
x_gate_damping = "data/amplitudes_single_qubit-X_mix_damping_prob-0.1.csv"
x_gate_damping_quest_probs = "data/probabilities_single_qubit-X_mix_damping_prob-0.1.csv"
damping_quest_probs = CSV.read(x_gate_damping_quest_probs,DataFrame,header=true)
prob = 0.1
amp_mat_damping = load_amp_vector_density_mat(x_gate_damping)
prob_vec_damping = get_prob_vec(amp_mat_damping)


# Simulated in Julia
ρ = Complex.([1,0])*transpose(Complex.([1,0]))
X = [0 1;1 0]
ρ̂ = @chain ρ begin
    update_ρ(_,X)
    apply_mix_damping(_,prob)
end

prob_vec_damping_sim = get_prob_vec(ρ̂)

@chain damping_quest_probs begin
    @rename :QuEST = $1
    @transform :JuliaQuEST = prob_vec_damping
    @transform :Julia = prob_vec_damping_sim
end



