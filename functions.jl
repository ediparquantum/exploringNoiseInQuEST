"""
    Compute rough order of memory needed to simulate N qubits
"""
function total_gq_needed_N_qubits(num_qubits,bits_per_float)
    bits_in_GB = 8e9
    num_float_per_qubit=4
    bits_per_qubit = bits_per_float*num_float_per_qubit
    bit_per_gb  = (bits_per_qubit/(bits_in_GB))
    total_gb_per_qubit = bit_per_gb*(2^num_qubits)
    return total_gb_per_qubit
end




function convert_double_vec_complex(path)
    df = @chain path begin
        CSV.read(_,DataFrame,header=true)
        @rename begin
            :imag = $2
        end 
        @rtransform :amp = Complex(:real,:imag) 
        @select :amp
    end
    return df.amp 
end

function reshape_into_density_matrix(amp_vec)
    dims = size(amp_vec,1) |> sqrt |> Int
    ρ = @chain amp_vec begin
        permutedims(_)
        reshape(_,dims,dims) 
    end
    return ρ 
end

function load_amp_vector_density_mat(path)
    return (reshape_into_density_matrix ∘ 
    convert_double_vec_complex)(path)
end

function get_prob_per_qubit(amp_mat)
    amp_probs = abs2.(amp_mat) |> normalize |> diag 
    @assert sum(amp_probs) == 1 "Total probability is not 1. Fix"
    return amp_probs
end


function update_ρ(ρ,gate)
    return gate*ρ*gate'
end


function get_prob_vec(ρ̂)
    
    if isdiag(ρ̂)
        prob_vec = Real.(ρ̂) |> diag
    else
        prob_vec = abs2.(ρ̂) |> diag
    end
    @assert abs(tr(ρ̂) - 1)<eps() "Total probability is not 1. Fix"
    return prob_vec
end


function compute_kraus_operators(prob)
    K₀ = [1 0;0 √(1-prob)]
    K₁ = [0 √(prob);0 0]
    return (K₀=K₀,K₁=K₁)
end

function apply_mix_damping(ρ,prob)
    kraus = compute_kraus_operators(prob)
    ρ_updated =  update_ρ(ρ,kraus[:K₀]) .+ update_ρ(ρ,kraus[:K₁])
    return ρ_updated
end




function enumerate_qubit_label_arrangement(N::Int)
    [bitstring((i))[end-N+1:end] for i ∈ 0:2^N-1]
end


function list_qubit_axis_arrangement(N::Int)
    arr = enumerate_qubit_label_arrangement(N)
    pattern = ["|"*i*"⟩" for i in arr]
    return pattern
end


function get_noise_prob(text)
    pattern = r"\d+(\.\d+)?"
    m = match(pattern, text)
    if isnothing(m)
        return "No noise"
    else
        return m.match
    end
end



"""
    FROM CHAT GPT - CHECK EXAMPLES WORK
    plot_probability_barplot(x_axis, probs, enum_comp_basis, noise_prob, basename)

Constructs a bar plot of circuit probabilities.

# Arguments
- `x_axis::Vector`: x-axis values for the bar plot.
- `probs::Vector`: probability values corresponding to the x-axis values.
- `enum_comp_basis::Vector`: labels for the x-axis ticks.
- `noise_prob::Float64`: noise probability used in the circuit.
- `basename::String`: name of the circuit.

# Returns
- `f::Figure`: the resulting plot.

# Details
- The function constructs a bar plot of circuit probabilities using the given input data.
- The `basename` is used to determine the legend label and gate name for the plot.
- If `basename` contains the string "damping", the legend label will include the noise probability.
- The function extracts the gate name from the `basename` using regular expression matching.
- If no match is found for "PauliX" or "HGate" in the `basename`, an empty string is used as the gate name.

# Example
```julia
x_axis = [1, 2, 3]
probs = [0.3, 0.4, 0.6]
enum_comp_basis = ["A", "B", "C"]
noise_prob = 0.1
basename = "prob-0.10.csv"

plot_probability_barplot(x_axis, probs, enum_comp_basis, noise_prob, basename)
"""
function plot_probability_barplot(x_axis,probs,enum_comp_basis,noise_prob,basename)
    
    if contains(basename,r"damping")
        leg_lab = "Circuit probabilities\nDamping, p=$(noise_prob)"
    else
        leg_lab = "Circuit probabilities"
    end

    mat = match(r"PauliX|HGate|HCZGate",basename)
    

    if isnothing(mat)
        @warn "No match to PauliX or HGate, update search criteria, setting gate to empty string"
        gate = ""
    elseif mat.match == "PauliX"
        gate = "Pauli X"
    elseif mat.match == "HGate"
        gate = "Hadamard"
    elseif mat.match == "HCZGate"
        gate = "Hadamard Each + CZ"
    else
        @warn "No match to PauliX or HGate, update search criteria, setting gate to empty string"
        gate = ""
    end

    f=Figure()
    ax=Axis(f[1,1],xlabel=L"|\psi\rangle",ylabel = L"P(|\psi\rangle)",xticks = (x_axis, enum_comp_basis))
    barplot!(ax,x_axis,probs,width=0.5,label=gate)
    ylims!(0,1)
    axislegend(leg_lab,position=:lt)
    return f
end
