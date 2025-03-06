"""
    AbstractAdductIon{S, T} <: AbstractChemical

Abstract type for adduct ions with core chemical type `S` and adduct type `T`.
"""
abstract type AbstractAdductIon{S, T} <: AbstractChemical end

"""
    AdductIon{S <: AbstractChemical, T <: AbstractAdduct} <: AbstractAdductIon{S, T}

Adduct ions forming in mass spectrometry.

# Fields
* `core`: the original chemical undergoing ionization. 
* `adduct`: the adduct of ion.

# Constructors
* `AdductIon(name::AbstractString, a::AbstractString; kwargs...)` -> `AdductIon(parse_chemical(name; kwargs...), parse_adduct(a))`
* `AdductIon(::Type{S}, name::AbstractString, a::AbstractString; kwargs...)` -> `AdductIon(parse_chemical(S, name; kwargs...), parse_adduct(a))`
* `AdductIon(cc::AbstractChemical, a::AbstractString)` -> `AdductIon(cc, parse_adduct(a))`
* `AdductIon(name::AbstractString, a::AbstractAdduct; kwargs...)` -> `AdductIon(parse_chemical(name; kwargs...), a)`
* `AdductIon(::Type{S}, name::AbstractString, a::AbstractAdduct; kwargs...)` -> `AdductIon(parse_chemical(S, name; kwargs...), a)`
"""
struct AdductIon{S <: AbstractChemical, T <: AbstractAdduct} <: AbstractAdductIon{S, T}
    core::S
    adduct::T
end
# """
#     ISF{S <: AbstractChemical, T <: AbstractAdduct} <: AbstractAdductIon{S, T}

# AdductIon forming in the ion source (in-source fragmentation).
# """
# struct ISF{S <: AbstractChemical, T <: AbstractAdduct} <: AbstractAdductIon{S, T}
#     core::S
#     adduct::T
# end

AdductIon(name::AbstractString, a::AbstractString; kwargs...) = AdductIon(parse_chemical(name; kwargs...), parse_adduct(a))
AdductIon(::Type{S}, name::AbstractString, a::AbstractString; kwargs...) where S = AdductIon(parse_chemical(S, name; kwargs...), parse_adduct(a))
AdductIon(cc::AbstractChemical, a::AbstractString) = AdductIon(cc, parse_adduct(a))
AdductIon(name::AbstractString, a::AbstractAdduct; kwargs...) = AdductIon(parse_chemical(name; kwargs...), a)
AdductIon(::Type{S}, name::AbstractString, a::AbstractAdduct; kwargs...) where S = AdductIon(parse_chemical(S, name; kwargs...), a)