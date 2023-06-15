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


"""
    get_noise_prob(text)

The `get_noise_prob` function extracts the noise probability from the given `text` using a regular expression pattern. It searches for a decimal number with or without a fractional part and returns the matched string.

# Arguments
- `text::AbstractString`: The input text from which to extract the noise probability.

# Returns
- `noise_prob::Union{AbstractString, Nothing}`: The extracted noise probability if found, or `nothing` if not found.

# Example
```julia
text = "prob-0.1.csv"
noise_prob = get_noise_prob(text)
# Output: "0.1"

"""
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
function plot_probability_barplot(x_axis,probs,enum_comp_basis,noise_prob,gate_label,noise_label)
    
    if noise_label != "Noiseless"
        leg_lab = "Circuit probabilities\n$(noise_label), p=$(noise_prob)"
    else
        leg_lab = "Circuit probabilities"
    end


    f=Figure()
    ax=Axis(f[1,1],xlabel=L"|\psi\rangle",ylabel = L"P(|\psi\rangle)",xticks = (x_axis, enum_comp_basis))
    barplot!(ax,x_axis,probs,width=0.5,label=gate_label)
    ylims!(0,1)
    axislegend(leg_lab,position=:lt)
    return f
end





"""
    rjoin_or(text::Vector)

Converts an array of strings into a regular expression pattern separated by the OR operator ("|").

## Arguments
- `text::Vector`: An array of strings.

## Returns
A `Regex` object representing the regular expression pattern.

## Examples
```julia
data_sub_folders = ["damping", "no_noise"]
regex_expression = rjoin_or(data_sub_folders)
"""
function rjoin_or(text::Vector)
    return join(text, "|") |> Regex
end


"""
    whitespace_removal(text::String)

Remove all whitespace characters from the given text.

# Arguments
- `text::String`: The input text to remove whitespace from.

# Returns
- The input text with all whitespace characters removed.

# Examples
```julia
julia> whitespace_removal("  Hello    World  ")
"HelloWorld"

julia> whitespace_removal("   This   is   a   sentence.   ")
"Thisisasentence."

"""
function whitespace_removal(text::String)
    text = replace(text," "=>"")
    return text
end




"""
    make_camel_case(text::String)

Convert the input text to camel case format.

# Arguments
- `text::String`: The input text to convert to camel case.

# Returns
- The input text converted to camel case format.

# Examples
```julia
julia> make_camel_case("hello world")
"HelloWorld"

julia> make_camel_case("this is a sentence")
"ThisIsASentence"

"""
function make_camel_case(text::String)
    return text |> titlecase |> whitespace_removal
end



    


"""
    unmake_camel_case(text)

Converts a camel case string to a space-separated string.

# Arguments
- `text::String`: The input string in camel case format.

# Returns
- `result::String`: The converted string with spaces separating words.

# Examples
```julia
julia> unmake_camel_case("helloWorld")
"hello World"

julia> unmake_camel_case("thisIsATest")
"this Is A Test"

"""
function unmake_camel_case(text)
    isempty(text) && return ""
    
    result = string(text[1])  # Add the first character as-is
    ltext = 2:length(text)
    for i in ltext
        if isuppercase(text[i]) && !isuppercase(text[i-1])
            result *= " "
        end
        result *= string(text[i])
    end
    
    return result
end


"""
    get_label(data_sub_regex, basename, no_match_return)

Extracts a label from a string using a regular expression.

# Arguments
- `data_sub_regex::Regex`: The regular expression pattern to match against the `basename`.
- `basename::String`: The input string from which the label will be extracted.
- `no_match_return::String`: The value to return if no match is found.

# Returns
- `label::String`: The extracted label.

# Examples
```julia
julia> get_label(r"damping|no_noise", "damping_data.csv", "Unknown")
"Damping Data"

julia> get_label(r"damping|no_noise", "unknown_data.csv", "Unknown")
"Unknown"

"""
function get_label(data_sub_regex,basename,no_match_return) 
    mat = match(data_sub_regex,basename)
    
    if isnothing(mat)
        @warn "No match found returning user input $(no_match_return)"
        return no_match_return
    else
        return unmake_camel_case(mat.match) |> titlecase
    end
end




function get_plot_save_path(basename_data_path)
    if contains(basename_data_path,r"damping")
        folder = "damping"
        text = replace(basename_data_path,"csv"=>"png")
    else
        folder = "no_noise"
        text = replace(basename_data_path,"csv"=>"png")
    end
    path = string("/",folder,"/",text)
    return path
end







"""
    save_plot_verify(plot_str, plot)

The `save_plot_verify` function saves a plot to a file specified by `plot_str` and verifies the successful saving of the plot. It checks if the file exists using the `isfile` function and returns an appropriate message indicating the status.

# Arguments
- `plot_str::AbstractString`: The file path to save the plot.
- `plot::Plots.Plot`: The plot object to be saved.

# Returns
- `status::String`: The status message indicating the success or failure of the plot saving.

# Example
```julia
plot_str = "plot.png"
plot = plot(x, y)
status = save_plot_verify(plot_str, plot)
# Output: "Plot successfully saved to path"
"""
function save_plot_verify(plot_str,plot)
    save(plot_str,plot)
    if isfile(plot_str) 
        return "Plot successfully save to path"
    else
        @error "Plot did not save to file as the file was not found using `isfile`"
    end
end