"""
    chemicalname(chemical::T; kwargs...) 

Get name of `chemical`. It is equivalent to `getchemicalattr(chemical, :name; kwargs...)` except that it returns `string("::", T)` when name is not available.
"""
function chemicalname(cc::T; kwargs...) where {T <: AbstractChemical}
    result = getchemicalattr(cc, :name; kwargs...)
    isnothing(result) ? string("::", T) : result
end

"""
    chemicalformula(chemical::AbstractChemical; kwargs...)

Get formula of `chemical`. It is equivalent to `getchemicalattr(chemical, :formula; kwargs...)` except that it returns `""` when formula is not available.
"""
function chemicalformula(cc::AbstractChemical; kwargs...) 
    result = getchemicalattr(cc, :formula; kwargs...)
    isnothing(result) ? "" : result
end

"""
    chemicalelements(chemical::AbstractChemical; kwargs...)

Get elements of `chemical`. It is equivalent to `getchemicalattr(chemical, :elements; kwargs...)` except that it returns `Pair{String, Int}[]` when elements are not available.
"""
function chemicalelements(cc::AbstractChemical; kwargs...) 
    result = getchemicalattr(cc, :elements; kwargs...)
    isnothing(result) ? Pair{String, Int}[] : result
end

"""
    chemicalabbr(chemical::AbstractChemical; kwargs...)

Get abbreviation of `chemical`. It is equivalent to `getchemicalattr(chemical, :abbreviation; kwargs...)` except that it returns `chemicalname(chemical; kwargs...)` when abbreviation is not available.
"""
function chemicalabbr(cc::AbstractChemical; kwargs...) 
    result = getchemicalattr(cc, :abbreviation; kwargs...)
    isnothing(result) ? chemicalname(m; kwargs...) : result
end

"""
    chemicalsmiles(chemical::AbstractChemical; kwargs...)

Get SMILES of `chemical`. It is equivalent to `getchemicalattr(chemical, :SMILES; kwargs...)` except that it returns `""` when SMILES is not available.
"""
function chemicalsmiles(cc::AbstractChemical; kwargs...) 
    result = getchemicalattr(cc, :SMILES; kwargs...)
    isnothing(result) ? "" : result
end

"""
    ioncore(adduct_ion::AbstractAdductIon; kwargs...)

Core chemical of `adduct_ion`.
"""
ioncore(adduct_ion::AbstractAdductIon; kwargs...) = getchemicalattr(adduct_ion, :core; kwargs...)

"""
    ionadduct(adduct_ion::AbstractAdductIon; kwargs...)

Adduct of `adduct_ion`.
"""
ionadduct(adduct_ion::AbstractAdductIon; kwargs...) = getchemicalattr(adduct_ion, :adduct; kwargs...)

"""
    kmer(adduct_ion::AbstractAdductIon; kwargs...)

The number of core chemical. For instance, 2 for "[2M+H]+".
"""
kmer(adduct_ion::AbstractAdductIon; kwargs...) = getchemicalattr(adduct_ion, :kmer; kwargs...)

"""
    charge(chemical::AbstractChemical; kwargs...)

The charge of `chemical` (positive or negative). For instance, -2 for dianion, +3 for trication. The default value for cstion and anion are 1 and -1.
"""
charge(cc::AbstractChemical; kwargs...) = getchemicalattr(cc, :charge; kwargs...)

"""
    ncharge(chemical::AbstractChemical; kwargs...)

The number of charges of `chemical`. 
"""
ncharge(cc::AbstractChemical; kwargs...) = abs(charge(cc; kwargs...))

"""
    abundantchemical(chemical::AbstractChemical)

The most abundant chemical from a chemical (itself) or isobars. 
"""
abundantchemical(cc::AbstractChemical) = getchemicalattr(cc, :charge; kwargs...)

"""
    rt(chemical::AbstractChemical; kwargs...)

Get retention time of `chemical`. It is equivalent to `getchemicalattr(chemical, :rt; kwargs...)` except that it returns `NaN` when rt is not available.
"""
function rt(cc::AbstractChemical; kwargs...) 
    result = getchemicalattr(cc, :rt; kwargs...)
    isnothing(result) ? NaN : result 
end
