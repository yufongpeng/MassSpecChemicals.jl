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
* `chemicals`: a vector of chemicals.
* `abundnace`: the abundance of each chemicals.
"""
struct Isobars{T <: AbstractChemical} <: AbstractChemical
    chemicals::Vector{T}
    abundance::Vector{Float64}
end

"""
    Isotopomers{T <: AbstractChemical} <: AbstractChemical

Chemicals differed from isotopic replacement location.

# Fields 
* `parent`: shared chemical structure of isotopomers prior to isotopic replacement. 
* `isotopes`: `Vector{Pair{String, Int}}`. Isotopes-number pairs of isotopic replacement.

# Constructors
* `Isotopomers(parent::AbstractChemical, fullformula::String)`
* `Isotopomers(parent::AbstractChemical, fullelements::Dictionary)`
"""
struct Isotopomers{T <: AbstractChemical} <: AbstractChemical
    parent::T 
    isotopes::Vector{Pair{String, Int}}
end

function Isotopomers(parent::AbstractChemical, fullformula::String)
    Isotopomers(parent, unique_elements(chemicalelements(fullformula)))
end

function Isotopomers(parent::AbstractChemical, fullelements::Dictionary)
    dp = unique_elements(chemicalelements(parent))
    dr = deepcopy(fullelements)
    for k in keys(fullelements)
        haskey(ELEMENTS[:ISOTOPES], k) && (delete!(dr, k); continue)
        dr[k] -= get(dp, k, 0) 
    end
    Isotopomers(parent, [k => v for (k, v) in pairs(dr)])
end

"""
    ChemicalLoss{T <: AbstractChemical} <: AbstractChemical

Chemical loss from a precursor. This product is not detected in MS; the other part of precursor is detected instead.
"""
struct ChemicalLoss{T <: AbstractChemical} <: AbstractChemical
    chemical::T 
end

"""
    ChemicalPair{T <: AbstractChemical, S <: AbstractChemical} <: AbstractChemical

A pair of precursor and product in MS/MS. Products can also be a `ChemicalLoss`.
"""
struct ChemicalPair{T <: AbstractChemical, S <: AbstractChemical} <: AbstractChemical
    precursor::T
    product::S
end