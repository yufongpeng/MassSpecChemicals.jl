"""
    AbstractIon{S, T} <: AbstractChemical

Abstract type for Ions.
"""
abstract type AbstractIon{S, T} <: AbstractChemical end

"""
    Ion{S <: AbstractChemical, T <: AbstractAdduct} <: AbstractIon{S, T}

Ion forming in the collisional cell.
"""
struct Ion{S <: AbstractChemical, T <: AbstractAdduct} <: AbstractIon{S, T}
    core::S
    adduct::T
end
"""
    ISF{S <: AbstractChemical, T <: AbstractAdduct} <: AbstractIon{S, T}

Ion forming in the ion source (in-source fragmentation).
"""
struct ISF{S <: AbstractChemical, T <: AbstractAdduct} <: AbstractIon{S, T}
    core::S
    adduct::T
end
"""
    IonCluster{T <: AbstractIon}

Ions with similar m/z.
"""
struct IonCluster{T <: AbstractIon}
    ions::Vector{T}
    abundance::Vector{Float64}
end

Ion(m::AbstractString, a::AbstractString; kwargs...) = Ion(Chemical("", m, collect(kwargs)), parse_adduct(a))
Ion(m::AbstractChemical, a::AbstractString) = Ion(m, parse_adduct(a))
Ion(m::AbstractString, a::AbstractAdduct; kwargs...) = Ion(Chemical("", m, collect(kwargs)), a)