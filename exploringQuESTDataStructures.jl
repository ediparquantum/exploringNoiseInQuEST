using Pkg
Pkg.activate(".")
using LinearAlgebra
using CSV
using DelimitedFiles
using DataFrames
using DataFramesMeta
using Chain
using CairoMakie
using BenchmarkTools
using Formatting
attributes = Theme(fontsize = 25,resolution=(800,800),Axis = (; aspect = 1),Legend = (;framevisible=false))
set_theme!(attributes)


// # Memory needed
    """
        Plot to show how much memory is needed in simulating qubits
    """
    num_qubits = 1:40
    bits_per_float=64
    gb = total_gq_needed_N_qubits(num_qubits,bits_per_float)
    f=Figure()
    ax=Axis(f[1,1],xlabel="Number of Qubits",ylabel="Gigabytes")
    scatter!(ax,qub,gb,markersize=10,label="GB to simulate")
    axislegend()
    f
    save("figs/memory_simulate_qubits.png",f)

//


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
