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
attributes = Theme(fontsize = 25,resolution=(800,800),Axis = (; aspect = 1),Legend = (;framevisible=false))
set_theme!(attributes)


state_vector_path = "data/density_matrix.csv"
h_state_vector_path = "data/density_matrix1.csv"
sv_df = CSV.read(state_vector_path,DataFrame,header=true)
h_sv_df = CSV.read(h_state_vector_path,DataFrame,header=true)



function convert_double_vec_complexDLM(path)
vec = @chain path begin
    readdlm(_,',') 
    _[2:end,:]
    Float64.(_)
    [Complex(i[1],i[2]) for i in eachrow(_)]
end
return vec
end

function convert_double_vec_complexDLMReinterpret(path)
    vec = @chain path begin
        readdlm(_,',') 
        _[2:end,:]
        Float64.(_)
        reinterpret(ComplexF64, _)
    end
    return vec
    end
    
    function convert_double_vec_complexDLMReinterpretNoChain(path)
        vec = readdlm(path,',')[2:end,:]
        vec = Float64.(vec)
        vec = reinterpret(ComplexF64, vec)
        
        return vec
    end


function convert_double_vec_complexCSV(path)
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
    return reshape(amp_vec,dims,dims)'
end

function load_amp_vector_density_mat(path)
    return (reshape_into_density_matrix ∘ 
        convert_double_vec_complexCSV)(path)
end

function get_prob_per_qubit(amp_mat)
    amp_probs = abs2.(amp_mat) |> normalize |> diag 
    @assert sum(amp_probs) == 1 "Total probability is not 1. Fixe"
    return amp_probs
end



amp_mat = load_amp_vector_density_mat(state_vector_path)
qub_probs = get_prob_per_qubit(amp_mat)


memGBqub(qub)=(256/(8e9))*(2^qub)

qub = 1:40
gb = memGBqub.(qub)
f=Figure()
ax=Axis(f[1,1],xlabel="Number of Qubits",ylabel="Gigabytes")
scatter!(ax,qub,gb,markersize=10,label="GB to simulate")
axislegend()
f
save("figs/memory_simulate_qubits.png",f)