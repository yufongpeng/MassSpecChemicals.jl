"""
    Chemical <: AbstractChemical

Unstructured chemical type with its name, formula, and additional information.
"""
struct Chemical <: AbstractChemical
    name::String
    formula::String
    info::Vector{Pair{Symbol, Any}}
end

Chemical(name::AbstractString, formula::AbstractString; kwargs...) = Chemical(name, formula, collect(kwargs))