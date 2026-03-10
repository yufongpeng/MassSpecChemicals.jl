"""
    AbstractAdductIon{S, T} <: AbstractChemical

Abstract type for adduct ions with core chemical type `S` and adduct type `T`.

# Special attributes
* `kmer -> Int`: number of core chemical "M" in adduct ion representation "[M+X]n+". 
* `ioncore -> S`: the core chemical undergoing ionization. 
* `ionadduct -> T`: the adduct formed during ionization. 
"""
abstract type AbstractAdductIon{S, T} <: AbstractChemical end

"""
    AdductIon{S <: AbstractChemical, T <: AbstractAdduct} <: AbstractAdductIon{S, T}

Adduct ions forming in mass spectrometry.

# Fields
* `core`: the core chemical undergoing ionization. 
* `adduct`: the adduct formed during ionization.

# Constructors
* `AdductIon(name::AbstractString, a::AbstractString; kwargs...)` -> `AdductIon(parse_chemical(name; kwargs...), parse_adduct(a))`
* `AdductIon(name::AbstractString, a::AbstractAdduct; kwargs...)` -> `AdductIon(parse_chemical(name; kwargs...), a)`
* `AdductIon(::Type{S}, name::AbstractString, a::AbstractString; kwargs...)` -> `AdductIon(parse_chemical(S, name; kwargs...), parse_adduct(a))`
* `AdductIon(::Type{S}, name::AbstractString, a::AbstractAdduct; kwargs...)` -> `AdductIon(parse_chemical(S, name; kwargs...), a)`
* `AdductIon(cc::AbstractChemical, a::AbstractString)` -> `AdductIon(cc, parse_adduct(a))`
"""
struct AdductIon{S <: AbstractChemical, T <: AbstractAdduct} <: AbstractAdductIon{S, T}
    core::S
    adduct::T
end

AdductIon(name::AbstractString, a::AbstractString; kwargs...) = AdductIon(parse_chemical(name; kwargs...), parse_adduct(a))
AdductIon(name::AbstractString, a::AbstractAdduct; kwargs...) = AdductIon(parse_chemical(name; kwargs...), a)
AdductIon(::Type{S}, name::AbstractString, a::AbstractString; kwargs...) where S = AdductIon(parse_chemical(S, name; kwargs...), parse_adduct(a))
AdductIon(::Type{S}, name::AbstractString, a::AbstractAdduct; kwargs...) where S = AdductIon(parse_chemical(S, name; kwargs...), a)
AdductIon(cc::AbstractChemical, a::AbstractString) = AdductIon(cc, parse_adduct(a))