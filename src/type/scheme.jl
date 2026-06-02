
"""
"""
abstract type AbstractElementalScheme <: AbstractScheme end
abstract type AbstractStructuralScheme <: AbstractScheme end
abstract type AbstractCompleteScheme{T, S} <: AbstractScheme end

struct RandomProductScheme <: AbstractStructuralScheme end
  
struct StructuralElementalScheme{T, S} <: AbstractCompleteScheme{T, S}
    structuralscheme::T 
    elementalscheme::S
end

"""

    ElementalScheme{Bool, T<:AbstractChemical} <: AbstractElementalScheme

Single scheme involving a chemical. The elements are fixed, and can be replaced by minor isotopes. 

* `ElementalScheme{false}`: chemical loss from a precursor. This product is not detected in MS; the other part of precursor is detected instead.
* `ElementalScheme{true}`: chemical gain to a precursor. This product is not detected in MS; the merged chemical is detected instead.

# Fields 
* `chemical::T`: chemical involved in scheme.

Single isotopomer can be set by using `Isotopomers` as field `chemical`. Elemental scheme can be redirected to correponding isotopic labeled scheme in `AdductIon` by dispatching on core chemical and existed schema, or looking up the `property` for generic `Chemical`.
"""
struct ElementalScheme{Bool, T<:AbstractChemical} <: AbstractElementalScheme
    chemical::T
    function ElementalScheme(gain::Bool, x::T) where {T <: AbstractChemical}
        new{gain, T}(x)
    end
end

"""
    ChemicalGain(chemical)

Chemical gain of `chemical`, i.e. `ElementalScheme(true, chemical)`.
"""
ChemicalGain(x) = ElementalScheme(true, x)

"""
    ChemicalLoss(chemical)

Chemical loss of `chemical`, i.e. `ElementalScheme(false, chemical)`.
"""
ChemicalLoss(x) = ElementalScheme(false, x)

"""
    ChemicalSchema{T<:AbstractScheme} <: AbstractScheme

Mutiple chemical schema.

# Fields
* `schema::Dictionary{T, Int}`: scheme => number dictionary
"""
struct ChemicalSchema{T<:AbstractScheme} <: AbstractScheme
    schema::Dictionary{T, Int}
end 

"""
    IsotopomerizedSchema{T<:AbstractScheme} <: AbstractScheme

Mutiple chemical schema with delocalized isotopic replacements.

# Fields
* `schema::ChemicalSchema{T}`: parent scheme
* `isotopes::Vector{Pair{String, Int}}`: delocalized isotopic replacements
"""
struct IsotopomerizedSchema{T<:AbstractScheme} <: AbstractScheme 
    schema::ChemicalSchema{T}
    isotopes::Vector{Pair{String, Int}}
end

schemetype(::IsotopomerizedSchema{T}) where T = T 
schemetype(::ChemicalSchema{T}) where T = T 
schemetype(::T) where T = T 

function ChemicalSchema(scheme::T, schema...) where {T <: AbstractScheme} 
    C = promote_type(T, schemetype.(schema)...)
    cs = Dictionary{C, Int}([scheme], [1])
    for s in schema
        push_scheme!(cs, s)
    end
    ChemicalSchema(cs)
end

function ChemicalSchema(schema::AbstractVector{T}) where {T <: AbstractScheme}
    cs = Dictionary{T, Int}([first(schema)], [1])
    length(schema) < 2 && return ChemicalSchema(cs)
    for s in @view schema[2:end]
        push_scheme!(cs, s)
    end
    ChemicalSchema(cs)
end

function ChemicalSchema(scheme::ChemicalSchema{T}, schema...) where {T <: AbstractScheme}
    C = promote_type(T, schemetype.(schema)...)
    if C == T
        cs = copy(scheme.schema)
    else
        cs = convert(Dictionary{C, Int}, copy(scheme.schema))
    end
    for s in schema
        push_scheme!(cs, s)
    end
    ChemicalSchema(cs)
end

function push_scheme!(cs::Dictionary, scheme::ChemicalSchema)
    for (k, v) in pairs(scheme.schema)
        get!(cs, k, 0)
        cs[k] += v
    end
    cs
end

function push_scheme!(cs::Dictionary, scheme::AbstractScheme)
    get!(cs, scheme, 0)
    cs[scheme] += 1
    cs
end

"""
    const CompleteSchema = Union{<:AbstractCompleteScheme, <:ChemicalSchema{<:AbstractCompleteScheme}, <:IsotopomerizedSchema{<:AbstractCompleteScheme}}

Complete scheme (scheme containing both `structuralscheme` and `elementalscheme`)
"""
const CompleteSchema = Union{<:AbstractCompleteScheme, <:ChemicalSchema{<:AbstractCompleteScheme}, <:IsotopomerizedSchema{<:AbstractCompleteScheme}}

"""
    const StructuralSchema = Union{<:AbstractStructuralScheme, <:ChemicalSchema{<:AbstractStructuralScheme}, <:IsotopomerizedSchema{<:AbstractStructuralScheme}}

Stuctural schema.
"""
const StructuralSchema = Union{<:AbstractStructuralScheme, <:ChemicalSchema{<:AbstractStructuralScheme}, <:IsotopomerizedSchema{<:AbstractStructuralScheme}}

"""
    const ElementalSchema = Union{<:AbstractElementalScheme, <:ChemicalSchema{<:AbstractElementalScheme}, <:IsotopomerizedSchema{<:AbstractElementalScheme}}

Elemental schema.
"""
const ElementalSchema = Union{<:AbstractElementalScheme, <:ChemicalSchema{<:AbstractElementalScheme}, <:IsotopomerizedSchema{<:AbstractElementalScheme}}
