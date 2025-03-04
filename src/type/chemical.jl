"""
    Chemical <: AbstractChemical

Unstructured chemical type with its name, formula, and additional attributes.

# Fields 
* `name`: a unique chemical name (`String`).
* `formula`: chemical formula (`String`).
* `attr`: additional attributes (`Vector{Pair{Symbol, Any}}`); the pairs repressent attribute names and values.

# Constructors
* `Chemical(name::AbstractString, formula::AbstractString; kwargs_as_attr_pairs...)`
"""
struct Chemical <: AbstractChemical
    name::String
    formula::String
    attr::Vector{Pair{Symbol, Any}}
end

Chemical(name::AbstractString, formula::AbstractString; kwargs...) = Chemical(name, formula, collect(kwargs))

"""
    Isobars{T <: AbstractChemical} <: AbstractChemical

Chemicals with similar m/z.

# Fields 
* `chemicals`: a vector of ions.
* `abundnace`: the abundance of each ions.
"""
struct Isobars{T <: AbstractChemical} <: AbstractChemical
    chemicals::Vector{T}
    abundance::Vector{Float64}
end