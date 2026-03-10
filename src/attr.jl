defaultname(::T) where {T <: AbstractChemical} = string("Chemical::", T)

"""
    chemicalname(chemical::AbstractChemical; verbose = true, kwargs...) -> String

The name of `chemical`. It should be unique for each chemical. 

# Methods
* `AbstractChemical`: property search
* `AbstractAdductIon`: generic function

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:name`
    * default: `defaultname(chemical)`

# Keyword arguments
* `verbose` determines whether includes all names or not for chemical species. If `verbose` is false, only the first (most abundant) chemical is included.
"""
chemicalname(chemical::AbstractChemical; kwargs...) = getchemicalproperty(chemical, :name, defaultname(chemical))

"""
    chemicalformula(chemical::AbstractChemical; kwargs...) -> String

The formula of chemical entity of `chemical`.

# Methods
* `AbstractChemical`: property search
* `AbstractAdductIon`: generic function

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:formula`
    * default: `""`
2. Function: `chemicalelements` -> `chemicalformula`
    * `kwargs` for `chemicalformula`. See documents of other methods.

# Keyword arguments
* `unique` determines whether combines the elements to become unique or not when constructing formula from attribute `chemicalelements`. It defaults to false for type `Chemical`.
"""
function chemicalformula(chemical::AbstractChemical; shallow = false, kwargs...) 
    result = getchemicalproperty(chemical, :formula, "")
    if isempty(result) && !shallow
        chemicalformula(chemicalelements(chemical; shallow = true); kwargs...)
    else 
        result 
    end
end

"""
    chemicalelements(chemical::AbstractChemical; kwargs...) -> Vector{Pair{String, Int}}

The elements of chemical entity of `chemical`.

# Methods
* `AbstractChemical`: property search
* `AbstractAdductIon`: generic function

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:elements`
    * default: `Pair{String, Int}[]`
2. Function: `chemicalformula` -> `chemicalelements`
    * `kwargs` for `chemicalelements`. See documents of other methods.
"""
function chemicalelements(chemical::AbstractChemical; shallow = false, kwargs...) 
    result = getchemicalproperty(chemical, :elements, Pair{String, Int}[])
    if isempty(result) && !shallow
        chemicalelements(chemicalformula(chemical; shallow = true); kwargs...)
    else 
        result 
    end
end

"""
    chemicalabbr(chemical::AbstractChemical; verbose = true, kwargs...) -> String

The abbreviation of `chemical`. 

# Methods
* `AbstractChemical`: property search
* `AbstractAdductIon`: generic function

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:abbreviation`
    * default: `chemicalname(chemical; kwargs...)`

# Keyword arguments
* `verbose` determines whether includes all abbreviations or not for chemical species. If `verbose` is false, only the first (most abundant) chemical is included.
"""
chemicalabbr(chemical::AbstractChemical; kwargs...) = getchemicalproperty(chemical, :abbreviation, chemicalname(chemical; kwargs...))

"""
    chemicalsmiles(chemical::AbstractChemical; kwargs...) -> String

The SMILES of chemical entity of `chemical`.

# Methods
* `AbstractChemical`: property search
* `AbstractAdductIon`: property search

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:SMILES`
    * default: `""`

For type `Isotopomers`, it returns SMILES of the parent.
For type `AbstractAdductIon`, it returns SMILES of the core chemical.
"""
chemicalsmiles(chemical::AbstractChemical; kwargs...) = getchemicalproperty(chemical, :SMILES, "")

"""
    ioncore(adduct_ion::AbstractAdductIon{S, T}; kwargs...) -> Union{S, Nothing}

The core chemical of `adduct_ion`. 

# Methods
* `AbstractAdductIon`: property search

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:core`
    * default: Error
"""
ioncore(adduct_ion::AbstractAdductIon{S, T}; kwargs...) where {S, T} = getchemicalproperty(adduct_ion, :core, nothing)::S

"""
    ionadduct(adduct_ion::AbstractAdductIon{S, T}; kwargs...) -> Union{T, Nothing}

The adduct of `adduct_ion`. 

# Methods
* `AbstractAdductIon`: property search

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:adduct`
    * default: Error
"""
ionadduct(adduct_ion::AbstractAdductIon{S, T}; kwargs...) where {S, T} = getchemicalproperty(adduct_ion, :adduct, nothing)::T

"""
    kmer(adduct_ion::AbstractAdductIon; kwargs...) -> Int

The number of core chemical. For instance, 2 for "[2M+H]+". 
"""
kmer(adduct_ion::AbstractAdductIon; kwargs...) = kmer(ionadduct(adduct_ion))

"""
    charge(chemical::AbstractChemical; kwargs...) -> Int

The charge state of `chemical`; positive for cation and negative for anion. For instance, -2 for dianion, +3 for trication. 

# Methods
* `AbstractChemical`: property search

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:charge`
    * default: `0`
"""
charge(chemical::AbstractChemical; kwargs...) = getchemicalproperty(chemical, :charge, 0)
charge(adduct_ion::AbstractAdductIon; kwargs...) = kmer(adduct_ion) * charge(ioncore(adduct_ion); kwargs...) + charge(ionadduct(adduct_ion))

"""
    ncharge(chemical::AbstractChemical; kwargs...) -> Int

The number of charges of `chemical`.
"""
ncharge(chemical::AbstractChemical; kwargs...) = abs(charge(chemical; kwargs...))

"""
    retentiontime(chemical::AbstractChemical; kwargs...) -> AbstractFloat

The retention time of `chemical`.

# Methods
* `AbstractChemical`: property search

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:retentiontime`
    * default: `NaN`
"""
retentiontime(chemical::AbstractChemical; kwargs...) = getchemicalproperty(chemical, :retentiontime, NaN)

"""
    chemicalentity(chemical::AbstractChemical) -> AbstractChemical

The single chemical entity from a chemical entity (i.e. itself), a chemical pair or chemical species.

By defaults, all chemical types are regarded as a chemical entity, except
* `Isobars`: `first(chemical.chemicals)`, i.e. the most abundant entity.
* `ChemicalLoss`: `chemical.chemical`.
* `ChemicalPair`: `chemicalentity(chemical.precursor)`, i.e. the very begining precursor.
"""
chemicalentity(chemical::AbstractChemical; kwargs...) = chemical

"""
    chemicalspecies(chemical::AbstractChemical) -> Vector{<: AbstractChemical}

The chemical species (represented by a vector) from a chemical entity (pair) or chemical species. 

By defaults, all chemical types are not chemical species, except
* `Isobars`: `chemicalchemicals`.
"""
chemicalspecies(chemical::AbstractChemical; kwargs...) = [chemical]

"""
    chemicalpair(chemical::AbstractChemical) -> Pair{<: AbstractChemical, <: AbstractChemical}

The chemical pair (represented by a pair). It returns `chemical => chemical` for `chemical` that is not paired.

By defaults, all chemical types are not chemical pairs, except
* `ChemicalPair`: `chemical.precursor => chemical.pair`.
"""
chemicalpair(chemical::AbstractChemical; kwargs...) = chemical => chemical

"""
    chemicalparent(chemical::AbstractChemical) -> AbstractChemical

The parent chemical without delocalized isotopes replacement.

By defaults, all chemical types are regarded as its own parent chemical, except
* `Isotopomers`: `chemical.parent`.
"""
chemicalparent(cc::AbstractChemical; kwargs...) = cc

"""
    isotopomersisotopes(chemical::AbstractChemical) -> Vector{Pair{String, Int}}

The delocalized isotopes replacement of isotopomers.

By defaults, all chemical types contain no replacement, except
* `Isotopomers`: `chemical.isotopes`.
"""
isotopomersisotopes(cc::AbstractChemical; kwargs...) = getchemicalproperty(cc, :isotopomersisotopes, Pair{String, Int}[])

# truly formed chemical
"""
    analyzedchemical(chemical::AbstractChemical) -> AbstractChemical

The single chemical entity that is directly analyzed in the very begining of instrumental analysis.

By defaults, all chemical types can be directly analyzed, except
* `ChemicalLoss`: throw error.
* `ChemicalPair`: `analyzedchemical(chemical.precursor)`, i.e. the very begining precursor.
"""
analyzedchemical(cc::AbstractChemical; kwargs...) = cc

"""
    detectedchemical(chemical::AbstractChemical) -> AbstractChemical

The single chemical entity that is directly detected in the very ending of instrumental analysis.

By defaults, all chemical types can be directly detected, except
* `ChemicalLoss`: precursor exclueded `chemical.chemical` part. It requires a keyword argument `precursor`.
* `ChemicalPair`: baically `detectedchemical(chemical.product)`, i.e. the very ending product, but takes `ChemicalLoss` in consideration.
"""
detectedchemical(cc::AbstractChemical; kwargs...) = cc

"""
    detectedisotopes(chemical::AbstractChemical) -> Vector{Pair{String, Int}}

The delocalized isotopes replacement of detected chemical.
"""
detectedisotopes(cc::AbstractChemical; kwargs...) = isotopomersisotopes(cc; kwargs...)

"""
    detectedcharge(chemical::AbstractChemical) -> Int

The charge state of detected chemical.
"""
detectedcharge(cc::AbstractChemical; kwargs...) = charge(cc; kwargs...)

"""
    detectedelements(chemical::AbstractChemical) -> Int

The elements of detected chemical.
"""
detectedelements(cc::AbstractChemical; kwargs...) = chemicalelements(cc; kwargs...)

# MS representation of chemical
"""
    inputchemical(chemical::AbstractChemical) -> AbstractChemical

The single chemical representation that is the input of the very begining of instrumental analysis. 

It is equivalent to `analyzedchemical`, except `ChemicalLoss` is accepeted.
"""
inputchemical(cc::AbstractChemical; kwargs...) = cc

"""
    outputchemical(chemical::AbstractChemical) -> AbstractChemical

The single chemical representation that is the output of the very ending of instrumental analysis.

It is equivalent to `detectedchemical`, except `ChemicalLoss` is kept unchanged.
"""
outputchemical(cc::AbstractChemical; kwargs...) = cc

# truly formed chemical in a tandem MS
"""
    analyzedprecursor(chemical::AbstractChemical) -> AbstractChemical

The single chemical entity that is directly analyzed in the nearest instrumental analysis.

By defaults, all chemical types can be directly analyzed, and for `ChemicalPair`, it is equivalent to `detectedchemical(chemical.precursor)`.
"""
analyzedprecursor(cc::AbstractChemical; kwargs...) = cc

"""
    detectedproduct(chemical::AbstractChemical) -> AbstractChemical

The single chemical entity that is directly detected in the nearest instrumental analysis.

By defaults, all chemical types can be directly detected, and for `ChemicalPair`, it is equivalent to `analyzedchemical(chemical.product)`.
"""
detectedproduct(cc::AbstractChemical; kwargs...) = cc

"""
    detectedproductisotopes(chemical::AbstractChemical) -> Vector{Pair{String, Int}}

The delocalized isotopes replacement of detected product.
"""
detectedproductisotopes(cc::AbstractChemical; kwargs...) = isotopomersisotopes(cc; kwargs...)

"""
    detectedproductcharge(chemical::AbstractChemical) -> Int

The charge state of detected product.
"""
detectedproductcharge(cc::AbstractChemical; kwargs...) = charge(cc; kwargs...)

"""
    detectedproductelements(chemical::AbstractChemical) -> Int

The elements of detected product.
"""
detectedproductelements(cc::AbstractChemical; kwargs...) = chemicalelements(cc; kwargs...)

# MS representation of chemical
"""
    inputprecursor(chemical::AbstractChemical) -> AbstractChemical

The single chemical representation that is the input of the nearest instrumental analysis. 

It is equivalent to `analyzedprecursor`, except `ChemicalLoss` is accepeted.
"""
inputprecursor(cc::AbstractChemical; kwargs...) = cc

"""
    outputproduct(chemical::AbstractChemical) -> AbstractChemical

The single chemical representation that is the output of the nearest instrumental analysis.

It is equivalent to `detectedproduct`, except `ChemicalLoss` is kept unchanged.
"""
outputproduct(cc::AbstractChemical; kwargs...) = cc

"""
    msstage(chemical::AbstractChemical) -> Int

The stage of MS the chemical has been through.
"""
msstage(cc::AbstractChemical; n = 0) = n + 1

"""
    chemicaltransitions(chemical::AbstractChemical) -> Vector{<: AbstractChemical}

The chemical entity(s) that are directly analyzed in each stage of instrumental analysis.

By defaults, all chemical types are directly analyzed in one stage of instrumental analysis, except
* `ChemicalPair`: a vector of analyzed chemicals.
"""
chemicaltransitions(cc::AbstractChemical) = [cc]
# list all propertys
# add new propertys
function chemicalpropertynames() end