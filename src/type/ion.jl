"""
    AbstractIon{S, T} <: AbstractChemical

Abstract type for Ions with core chemical type `S` and adduct type `T`.
"""
abstract type AbstractIon{S, T} <: AbstractChemical end

"""
    Ion{S <: AbstractChemical, T <: AbstractAdduct} <: AbstractIon{S, T}

Ions forming in mass spectrometry.

# Fields
* `core`: the chemical undergoing ionization. 
* `adduct`: the adduct of ion.

# Constructors
* `Ion(core::AbstractChemical, adduct::AbstractAdduct)`
* `Ion(formula::AbstractString, adduct::AbstractString; kwargs...)` -> `Ion(parse_chemical(formula; kwargs...), parse_adduct(adduct))`
* `Ion(core::AbstractChemical, adduct::AbstractString)` -> `Ion(core, parse_adduct(adduct))`
* `Ion(formula::AbstractString, adduct::AbstractAdduct; kwargs...)` -> `Ion(parse_chemical(formula; kwargs...), adduct)`

"""
struct Ion{S <: AbstractChemical, T <: AbstractAdduct} <: AbstractIon{S, T}
    core::S
    adduct::T
end
# """
#     ISF{S <: AbstractChemical, T <: AbstractAdduct} <: AbstractIon{S, T}

# Ion forming in the ion source (in-source fragmentation).
# """
# struct ISF{S <: AbstractChemical, T <: AbstractAdduct} <: AbstractIon{S, T}
#     core::S
#     adduct::T
# end

Ion(m::AbstractString, f::AbstractString, a::AbstractString; kwargs...) = Ion(parse_chemical(m, f; kwargs...), parse_adduct(a))
Ion(m::AbstractChemical, a::AbstractString) = Ion(m, parse_adduct(a))
Ion(m::AbstractString, f::AbstractString, a::AbstractAdduct; kwargs...) = Ion(parse_chemical(m, f; kwargs...), a)