module BioChemicals

using ..MassSpecChemicals
using ..MassSpecChemicals: AbstractChemical
using IterTools
export FunctionalGroup, LeavingGroup, Dehydrogen, Odehydrogen, Ndehydrogen, Didehydrogen, Odidehydrogen, Ndidehydrogen, Deamine, Dehydroxy, 
        Substituent, SubstitutedChemical, ChainedChemical, DehydratedChemical, IsotopiclabeledChemical,
        AbstractLinkageposition, Linkageposition, AbstractChainedChemical, 
        originalmolecule, leavinggroup, conjugation, includeSIL

abstract type FunctionalGroup{M, S} end
abstract type LeavingGroup end
struct Dehydrogen <: LeavingGroup end
struct Odehydrogen <: LeavingGroup end
struct Ndehydrogen <: LeavingGroup end
struct Didehydrogen <: LeavingGroup end
struct Odidehydrogen <: LeavingGroup end
struct Ndidehydrogen <: LeavingGroup end
struct Dehydroxy <: LeavingGroup end
struct Deamine <: LeavingGroup end
abstract type AbstractLinkageposition end
struct Linkageposition <: AbstractLinkageposition
    position::UInt8
end
lk(x) = Linkageposition(UInt8(x))

struct Substituent{M <: AbstractChemical, S <: LeavingGroup} <: FunctionalGroup{M, S}
    molecule::M
    position::AbstractLinkageposition
end
struct SubstitutedChemical{M <: AbstractChemical, S, T} <: AbstractChemical
    molecule::M
    substituted::S
    substituent::T
end 
#= 
S
Vector{Pair{LeavingGroup, UInt8}}: LeavingGroup => Number
Vector{Pair{AbstractLinkageposition, LeavingGroup}}: Position => LeavingGroup
T
Vector{Pair{FunctionalGroup, UInt8}}: FunctionalGroup => Number
Vector{Pair{AbstractLinkageposition, FunctionalGroup}}: Position => FunctionalGroup
=#
struct ChainedChemical{M, L} <: AbstractChemical
    molecule::M
    linkage::L
end
#=
L
Vector{Pair{LeavingGroup, LeavingGroup}}: LeavingGroup => LeavingGroup
Vector{Vector{Pair{LeavingGroup, LeavingGroup}}, Vector{Pair{AbstractLinkageposition, AbstractLinkageposition}}}: LeavingGroup => LeavingGroup Position([a/b]p for hex) => Position
=#
struct DehydratedChemical{M, L} <: AbstractChemical
    molecule::M
    linkage::L
end
#=
L
Nothing
Vector{Pair{AbstractLinkageposition, AbstractLinkageposition}}: Position([a/b]p for hex) => Position
=#
struct IsotopiclabeledChemical{M, I} <: AbstractChemical
    molecule::M
    isotopiclabel::I
end
Substituent{S}(m::M, p = UInt8(0)) where {M, S} = Substituent{M, S}(m, p)
originalmolecule(fg::Substituent) = fg.molecule
originalmolecule(sm::SubstitutedChemical) = sm.molecule
originalmolecule(::FunctionalGroup{M}) where M = M()
leavinggroup(::FunctionalGroup{M, S}) where {M, S} = S()
conjugation(m::AbstractChemical) = Substituent{Ndehydrogen}(m)
includeSIL(T) = Union{<: T, <: IsotopiclabeledChemical{<: T}}
const AbstractChainedChemical = [ChainedChemical, DehydratedChemical]

function makemolecule(m; sil)
    isnothing(sil) ? m : IsotopiclabeledChemical(m, s)
end
function makemolecule(::Type{T}, ms...; linkage = nothing, sil) where T
    sils = isnothing(sil) ? sil : distribute_sil(sil)
    molecules = ntuple(length(ms)) do i
        j = isnothing(sils) ? nothing : findfirst(x -> i == first(x), sils)
        isnothing(j) ? ms[i] : 
        in(typeof(ms[i]), AbstractChainedChemical) ? makemolecule(typeof(ms[i]), ms[i].chain...; linkage = ms[i].linkage, sil = last(sils[j])) : 
        IsotopiclabeledChemical(ms[i], last(sils[j]))
    end
    T(molecules, !isnothing(linkage) ? linkage : T == ChainedChemical ? [makelinkage(a, b) for (a, b) in IterTools.partition(molecules, 2, 1)] : nothing)
end
function distribute_sil(s)
    [parse(Int, a) => b for (a, b) in eachmatch(r"(\d)-(.*?[A-Z][a-z]*\d*)", s)]
end
# PC[1-D5] PC[3-1-(1,1,2,2)D4] PC[3-2-D9]
makelinkage(a, b) = Dehydroxy() => Dehydrogen()

include("basiccompound.jl")
include("metabolite.jl")
include(joinpath("glycan", "glycan.jl"))
include("aminoacid.jl")
include("nucleotide.jl")
include(joinpath("lipid", "lipid.jl"))

end