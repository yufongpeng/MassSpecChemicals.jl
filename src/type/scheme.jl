"""
    AbstractElementalScheme <: AbstractScheme end

Abstract scheme contaning exact elements.
"""
abstract type AbstractElementalScheme <: AbstractScheme end
"""
    AbstractStructuralScheme <: AbstractScheme end

Abstract scheme contaning only structures.  
"""
abstract type AbstractStructuralScheme <: AbstractScheme end
"""
    AbstractCompleteScheme{T, S} <: AbstractScheme end

Abstract scheme contaning both elements and structures. 
"""
abstract type AbstractCompleteScheme{T, S} <: AbstractScheme end
"""
    StructuralChemicalScheme <: AbstractStructuralScheme end

Abstract scheme contaning structures generating chemical entity. 
"""
abstract type StructuralChemicalScheme <: AbstractStructuralScheme end

struct RandomProductScheme <: AbstractStructuralScheme end

"""
    StructuralElementalScheme{T, S} <: AbstractCompleteScheme{T, S} end

Default [`AbstractCompleteScheme`](@ref).

# Fields
* `structuralscheme::T`.
* `elementalscheme::S`.
"""
struct StructuralElementalScheme{T, S} <: AbstractCompleteScheme{T, S}
    structuralscheme::T 
    elementalscheme::S
end

"""

    ElementalScheme{Bool, T<:AbstractChemical} <: AbstractElementalScheme
    ElementalScheme(gain::Bool, chemical::T) = ElementalScheme{gain, T}(chemical)
    
Single scheme involving a chemical. The elements are fixed, and can be replaced by minor isotopes. 
* `ElementalScheme{false}`: chemical loss from a precursor. This product is not detected in MS; the other part of precursor is detected instead.
* `ElementalScheme{true}`: chemical gain to a precursor. This product is not detected in MS; the merged chemical is detected instead.

# Fields 
* `chemical::T`: chemical involved in scheme.

Single isotopomer can be set by using [`Isotopomers`](@ref) as field `chemical`. 
Elemental scheme can be redirected to the corresponding isotopic labeled scheme in [`AdductIon`](@ref) by dispatching on core chemical and existing scheme, 
or looking up the `property` for generic [`Chemical`](@ref).
"""
struct ElementalScheme{Bool, T<:AbstractChemical} <: AbstractElementalScheme
    chemical::T
    function ElementalScheme(gain::Bool, x::T) where {T<:AbstractChemical}
        new{gain, T}(x)
    end
end

"""
    ChemicalGain(chemical)

Chemical gain of `chemical`, i.e. `ElementalScheme(true, chemical)`. See [`ElementalScheme`](@ref).
"""
ChemicalGain(x) = ElementalScheme(true, x)

"""
    ChemicalLoss(chemical)

Chemical loss of `chemical`, i.e. `ElementalScheme(false, chemical)`. See [`ElementalScheme`](@ref).
"""
ChemicalLoss(x) = ElementalScheme(false, x)

"""
    ChemicalSchemes{T<:AbstractScheme} <: AbstractScheme

Mutiple chemical schemes.

# Fields
* `schemes::Vector{Pair{T, Int}}`: scheme => number vector. The number v is the times of the scheme in the chemical schemes. 
"""
struct ChemicalSchemes{T<:AbstractScheme} <: AbstractScheme
    schemes::Vector{T}
    number::Vector{Int}
end 

"""
    IsotopomerizedSchemes{T<:ChemicalSchemes} <: AbstractScheme

Mutiple chemical schemes with delocalized isotopic replacements.

# Fields
* `schemes::T`: parent scheme.
* `isotopes::ElementsVector`: delocalized isotopic replacements.
"""
struct IsotopomerizedSchemes{T<:ChemicalSchemes} <: AbstractScheme 
    parent::T
    isotopes::ElementsVector
end

function IsotopomerizedSchemes(chemical::AbstractScheme, fullformula::String)
    IsotopomerizedSchemes(chemicalparent(chemical), dictionary_elements(chemicalelements(fullformula)))
end

function IsotopomerizedSchemes(chemical::AbstractScheme, fullelements::Vector{Pair{String, Int}})
    IsotopomerizedSchemes(chemicalparent(chemical), dictionary_elements(fullelements))
end

function IsotopomerizedSchemes(chemical::AbstractScheme, fullelements::Dict)
    parent = chemicalparent(chemical)
    dp = dictionary_elements(chemicalelements(parent))
    ev = ElementsVector(collect(keys(fullelements)), collect(values(fullelements)))
    del = Int[]
    for (i, (k, n)) in enumerate(ev)
        iselement(k) && (push!(del, i); continue)
        ev.numbers[i] = n - get(dp, k, 0) 
    end
    deleteat!(ev.elements, del)
    deleteat!(ev.numbers, del)
    IsotopomerizedSchemes(parent, ev)
end

"""
    Groupedisotopomerizedschemes{T<:AbstractScheme, N} <: AbstractScheme

Isotopomerized schemes grouped by mass-shift index.

# Fields 
* `parent::T`: shared chemical scheme prior to isotopic replacement. 
* `index::Int`: mass-shift index.
* `isotope::String`: reference isotope for [`mass_shift_index`](@ref).
* `isotopes::Vector{ElementsVector}`: Isotopes-number pairs of isotopic replacements of each isotopomers.
* `abundance::Vector{N}`: abundance of each isotopomers.
"""
struct Groupedisotopomerizedschemes{T<:AbstractScheme, N} <: AbstractScheme
    parent::T 
    index::Int
    isotope::String
    isotopes::Vector{ElementsVector}
    abundance::Vector{N}
    function Groupedisotopomerizedschemes(parent::T, index::Int, isotope::String, isotopes::Vector{ElementsVector}, abundance::Vector{N}) where {T, N}
        id = sortperm(abundance)
        new{T, N}(parent, index, isotope, isotopes[id], abundance[id])
    end
end

schemetype(::ChemicalSchemes{T}) where T = T 
schemetype(::T) where T = T 

function ChemicalSchemes(scheme::T, schemes...) where {T<:AbstractScheme} 
    C = promote_type(T, schemetype.(schemes)...)
    cs = C[scheme]
    cn = Int[1]
    for s in schemes
        push_scheme!(cs, cn, s)
    end
    ChemicalSchemes(cs, cn)
end

function ChemicalSchemes(schemes::AbstractVector{T}) where {T<:AbstractScheme}
    cs = T[first(schemes)]
    cn = Int[1]
    length(schemes) < 2 && return ChemicalSchemes(cs, cn)
    for s in @view schemes[2:end]
        push_scheme!(cs, cn, s)
    end
    ChemicalSchemes(cs, cn)
end

function ChemicalSchemes(scheme::ChemicalSchemes{T}, schemes...) where {T<:AbstractScheme}
    C = promote_type(T, schemetype.(schemes)...)
    if C == T
        cs = copy(scheme.schemes)
    else
        cs = convert(Vector{C}, copy(scheme.schemes))
    end
    cn = copy(scheme.number)
    for s in schemes
        push_scheme!(cs, cn, s)
    end
    ChemicalSchemes(cs, cn)
end

function push_scheme!(cs::Vector, cn::Vector, scheme::ChemicalSchemes)
    for (k, v) in zip(scheme.schemes, scheme.number)
        i = findfirst(==(k), cs)
        if i !== nothing
            cn[i] += v
        else
            push!(cs, k)
            push!(cn, v)
        end
    end
    cs
end

function push_scheme!(cs::Vector, cn::Vector, scheme::AbstractScheme)
    i = findfirst(==(scheme), cs)
    if i !== nothing
        cn[i] += 1
    else
        push!(cs, scheme)
        push!(cn, 1)
    end
    cs
end

"""
    const CompleteSchemes = Union{<:AbstractCompleteScheme, <:ChemicalSchemes{<:AbstractCompleteScheme}, <:IsotopomerizedSchemes{<:ChemicalSchemes{<:AbstractCompleteScheme}}, <:Groupedisotopomerizedschemes{<:ChemicalSchemes{<:AbstractCompleteScheme}}}

Complete scheme (scheme containing both [`structuralscheme`](@ref) and [`elementalscheme`](@ref)).
"""
const CompleteSchemes = Union{<:AbstractCompleteScheme, <:ChemicalSchemes{<:AbstractCompleteScheme}, <:IsotopomerizedSchemes{<:ChemicalSchemes{<:AbstractCompleteScheme}}, <:Groupedisotopomerizedschemes{<:ChemicalSchemes{<:AbstractCompleteScheme}}}

"""
    const StructuralSchemes = Union{<:AbstractStructuralScheme, <:ChemicalSchemes{<:AbstractStructuralScheme}, <:IsotopomerizedSchemes{<:ChemicalSchemes{<:AbstractStructuralScheme}}, <:Groupedisotopomerizedschemes{<:ChemicalSchemes{<:AbstractStructuralScheme}}}

Stuctural schemes.
"""
const StructuralSchemes = Union{<:AbstractStructuralScheme, <:ChemicalSchemes{<:AbstractStructuralScheme}, <:IsotopomerizedSchemes{<:ChemicalSchemes{<:AbstractStructuralScheme}}, <:Groupedisotopomerizedschemes{<:ChemicalSchemes{<:AbstractStructuralScheme}}}

"""
    const ElementalSchemes = Union{<:AbstractElementalScheme, <:ChemicalSchemes{<:AbstractElementalScheme}, <:IsotopomerizedSchemes{<:ChemicalSchemes{<:AbstractElementalScheme}}, <:Groupedisotopomerizedschemes{<:ChemicalSchemes{<:AbstractElementalScheme}}}

Elemental schemes.
"""
const ElementalSchemes = Union{<:AbstractElementalScheme, <:ChemicalSchemes{<:AbstractElementalScheme}, <:IsotopomerizedSchemes{<:ChemicalSchemes{<:AbstractElementalScheme}}, <:Groupedisotopomerizedschemes{<:ChemicalSchemes{<:AbstractElementalScheme}}}

@deprecate ChemicalSchema ChemicalSchemes
@deprecate IsotopomerizedSchema IsotopomerizedSchemes
@deprecate Groupedisotopomerizedschemes Groupedisotopomerizedschemes
@deprecate CompleteSchema CompleteSchemes
@deprecate StructuralSchema StructuralSchemes
@deprecate ElementalSchema ElementalSchemes

"""
    const CompleteSchemeChemical = AbstractCompleteScheme{T, <:AbstractChemical}

Complete scheme genrating chemical.
"""
const CompleteSchemeChemical = AbstractCompleteScheme{T, <:AbstractChemical} where T