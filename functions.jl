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