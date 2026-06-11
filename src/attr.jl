"""
    getchemicalproperty(chemical::AbstractChemicalsSchema, property::Symbol, default = nothing)

Get property from `chemical`. 

This function defaults to finds the property, and returns `default` If it is not available.

For type `Chemical` and properties other than `name`, `elements` and `formula`, it iterates through `chemical.property`. If no matched property name is found, it returns `default`.

For type `AbstractAdductIon`, it searches for properties of itself and then the properties of core chemical without specialized methods.
"""
getchemicalproperty(chemical::AbstractChemicalsSchema, property::Symbol, default = nothing) = hasproperty(chemical, property) ? getproperty(chemical, property) : default 
getchemicalproperty(adduct_ion::AbstractAdductIon, property::Symbol, default) = hasproperty(adduct_ion, property) ? getproperty(adduct_ion, property) : getchemicalproperty(ioncore(adduct_ion), property, default)
getchemicalproperty(chemical::AbstractChemicalWrapper, property::Symbol, default = nothing) = hasproperty(chemical, property) ? getproperty(chemical, property) : hasproperty(chemical.chemical, property) ? getproperty(chemical.chemical, property) : default 
function getchemicalproperty(chemical::Union{Chemical, FormulaChemical}, property::Symbol, default)
    hasproperty(chemical, property) && return getproperty(chemical, property)
    for (p, v) in chemical.property
        p == property && return v
    end
    return default
end

"""
    chemicalname(chemical::AbstractChemicalsSchema; verbose = true, n = 1, loss = false, bracket = true, kwargs...) -> String

The name of `chemical`. It should be unique for each object. 

# Generic Methods
* `AbstractChemicalsSchema`: property search
* `AbstractAdductIon`: name composed of name of core chemical and adduct
* `AbstractChemicalWrapper`: name of field `chemical`

# Specific Methods 
* Species/Transition Level

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:name`
    * default: `defaultname(chemical)`

# Keyword arguments
* `verbose` determines whether includes all names or not for chemical species. If `verbose` is false, only the first (most abundant) chemical is included.
* `n` determines number of chemical.
* `loss` determines whether the chemical is part of chemical loss, and sign flips are propagated into this chemical.
* `bracket` determines whether wrapping the output string with "[...]" and charge.
"""
function chemicalname(chemical::AbstractChemicalsSchema; n = 1, kwargs...) 
    nm = getchemicalproperty(chemical, :name, defaultname(chemical))
    n > 1 ? string(n, nm) : nm
end
chemicalname(chemical::AbstractChemicalWrapper; n = 1, kwargs...) = chemicalname(chemical.chemical; n, kwargs...)
function chemicalname(adduct_ion::AbstractAdductIon; loss = false, bracket = true, n = 1, kwargs...) 
    r = chemicalname(ioncore(adduct_ion); loss = bracket ? false : loss, bracket = false, n, kwargs...)
    a = chemicalabbr(ionadduct(adduct_ion); loss = bracket ? false : loss, bracket = false, n, kwargs...)
    n = ncore(adduct_ion)
    if occursin(" ", r)
        r = string("(", r, ")")
    end
    if n > 1 
        r = string(n, r)
    end
    bracket ? string("[", r, a, "]", charge_repr(charge(adduct_ion))) : string(r, a)
end

"""
    chemicalformula(chemical::AbstractChemicalsSchema; kwargs...) -> String

The formula of chemical entity of `chemical`.

# Generic Methods
* `AbstractChemical`: property search
* `AbstractAdductIon`: formuala combining core chemical and adduct
* `AbstractChemicalWrapper`: formula of field `chemical`
* `AbstractScheme`: property search; started with `"+"` for chemical gain and `"-"` for chemical loss

# Specific Methods 
* Entity Level

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:formula`
    * default: `""`
2. Function: `chemicalelements` -> `chemicalformula`

# Keyword arguments
* `delim` assigns the delimiter between each element.
* `unique` determines whether combines the elements to become unique or not when constructing formula from attribute `chemicalelements`. It defaults to false for type `Chemical` and `FormulaChemical`.
* `ischemical` determines whether the chemical is a chemical or a scheme. 
* `loss` determines whether the chemical is part of chemical loss, and signs are factored out from elements. 
"""
function chemicalformula(chemical::AbstractChemical; shallow = false, loss = false, ischemical = true, kwargs...) 
    result = getchemicalproperty(chemical, :formula, "")
    if isempty(result) && !shallow
        chemicalformula(chemicalelements(chemical; shallow = true); loss, ischemical, kwargs...)
    else 
        result 
    end
end

function chemicalformula(chemical::AbstractScheme; shallow = false, loss = true, ischemical = false, kwargs...) 
    result = getchemicalproperty(chemical, :formula, "")
    if isempty(result) && !shallow
        chemicalformula(chemicalelements(chemical; shallow = true); loss, ischemical, kwargs...)
    else 
        result 
    end
end

function chemicalformula(adduct_ion::AbstractAdductIon; kwargs...)
    el = copy(chemicalelements(ioncore(adduct_ion); kwargs...))
    nm = ncore(adduct_ion)
    if nm > 1
        for id in eachindex(el)
            el[id] = first(el[id]) => nm * last(el[id])
        end
    end
    chemicalformula(gain_elements!(dictionary_elements(Dictionary, el), chemicalelements(ionadduct(adduct_ion); kwargs..., loss = false)); kwargs...)
end

chemicalformula(chemical::AbstractChemicalWrapper; kwargs...) = chemicalformula(chemical.chemical; kwargs...)

"""
    chemicalelements(chemical::AbstractChemicalsSchema; kwargs...) -> Vector{Pair{String, Int}}

The elements of chemical entity of `chemical`.

# Generic Methods
* `AbstractChemical`: property search
* `AbstractAdductIon`: elements combining core chemical and adduct
* `AbstractChemicalWrapper`: chemical entity of field `chemical`
* `AbstractScheme`: property search; loss of elements

# Specific Methods 
* Entity Level

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:elements`
    * default: `Pair{String, Int}[]`
2. Function: `chemicalformula` -> `chemicalelements`

# Keyword Arguments 
* `loss` determines whether the chemical is part of chemical loss, and sign flips are propagated into elements.
"""
function chemicalelements(chemical::AbstractChemical; shallow = false, loss = false, kwargs...) 
    result = getchemicalproperty(chemical, :elements, Pair{String, Int}[])
    if isempty(result) && !shallow
        chemicalelements(chemicalformula(chemical; shallow = true); loss, kwargs...)
    else 
        result 
    end
end

function chemicalelements(chemical::AbstractScheme; shallow = false, loss = true, kwargs...) 
    result = getchemicalproperty(chemical, :elements, Pair{String, Int}[])
    if isempty(result) && !shallow
        chemicalelements(chemicalformula(chemical; shallow = true); loss, kwargs...)
    else 
        result 
    end
end

function chemicalelements(adduct_ion::AbstractAdductIon; kwargs...) 
    el = copy(chemicalelements(ioncore(adduct_ion); kwargs...))
    nm = ncore(adduct_ion)
    if nm > 1
        for id in eachindex(el)
            el[id] = first(el[id]) => nm * last(el[id])
        end
    end
    vcat(el, chemicalelements(ionadduct(adduct_ion); kwargs..., loss = false))
end

chemicalelements(chemical::AbstractChemicalWrapper; kwargs...) = chemicalelements(chemical.chemical; kwargs...)

"""
    chemicalabbr(chemical::AbstractChemicalsSchema; verbose = true, n = 1, loss = false, bracket = true, kwargs...) -> String

The abbreviation of `chemical`. 

# Target Chemical Level
* Species

# Generic Methods
* `AbstractChemicalsSchema`: property search
* `AbstractAdductIon`: abbreviation composed of abbreviation of core chemical and adduct
* `AbstractChemicalWrapper`: abbreviation of field `chemical`

# Specific Methods 
* Species/Transition Level

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:abbreviation`
    * default: `chemicalname(chemical; kwargs...)`

# Keyword arguments
* `verbose` determines whether includes all abbreviations or not for chemical species. If `verbose` is false, only the first (most abundant) chemical is included.
* `n` determines number of chemical.
* `loss` determines whether the chemical is part of chemical loss, and sign flips are propagated into this chemical.
* `bracket` determines whether wrapping the output string with "[...]" and charge.
"""
function chemicalabbr(chemical::AbstractChemicalsSchema; n = 1, kwargs...) 
    abbr = getchemicalproperty(chemical, :abbreviation, "")
    isempty(abbr) && return chemicalname(chemical; n, kwargs...)
    string(n > 1 ? n : "", abbr)
end
chemicalabbr(chemical::AbstractChemicalWrapper; n = 1, kwargs...) = chemicalabbr(chemical.chemical; n, kwargs...)
function chemicalabbr(adduct_ion::AbstractAdductIon; loss = false, bracket = true, n = 1, kwargs...) 
    r = chemicalabbr(ioncore(adduct_ion); loss = bracket ? false : loss, bracket = false, n, kwargs...)
    a = chemicalabbr(ionadduct(adduct_ion); loss = bracket ? false : loss, bracket = false, n, kwargs...)
    n = ncore(adduct_ion)
    if occursin(" ", r)
        r = string("(", r, ")")
    end
    if n > 1 
        r = string(n, r)
    end
    bracket ? string("[", r, a, "]", charge_repr(charge(adduct_ion))) : string(r, a)
end

"""
    chemicalsmiles(chemical::AbstractChemical; kwargs...) -> String

The SMILES of chemical entity of `chemical`.

# Generic Methods
* `AbstractChemical`: property search
* `AbstractAdductIon`: SMILES of the core chemical
* `AbstractChemicalWrapper`: SMILES of field `chemical`

# Specific Methods 
* Entity Level
* `Isotopomers`: SMILES of the parent chemical.
* `Groupedisotopomers`: SMILES of the most abundant parent chemical

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:SMILES`
    * default: `""`
"""
chemicalsmiles(chemical::AbstractChemical; kwargs...) = getchemicalproperty(chemical, :SMILES, "")
chemicalsmiles(adduct_ion::AbstractAdductIon; kwargs...) = chemicalsmiles(ioncore(adduct_ion); kwargs...) 
chemicalsmiles(chemical::AbstractChemicalWrapper; kwargs...) = chemicalsmiles(chemical.chemical; kwargs...)

"""
    ioncore(adduct_ion::AbstractAdductIon{S, T}; kwargs...) -> Union{S, Nothing}

The core chemical of `adduct_ion`. 

# Generic Methods
* `AbstractAdductIon`: property search

# Specific Methods 
* Entity Level

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:core`
    * default: Error
"""
ioncore(adduct_ion::AbstractAdductIon{S, T}; kwargs...) where {S, T} = getchemicalproperty(adduct_ion, :core, nothing)::S

"""
    ionadduct(adduct_ion::AbstractAdductIon{S, T}; kwargs...) -> Union{T, Nothing}

The adduct of `adduct_ion`. 

# Generic Methods
* `AbstractAdductIon`: property search

# Specific Methods 
* Entity Level

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:adduct`
    * default: Error
"""
ionadduct(adduct_ion::AbstractAdductIon{S, T}; kwargs...) where {S, T} = getchemicalproperty(adduct_ion, :adduct, nothing)::T

"""
    ncore(adduct_ion::AbstractAdductIon; kwargs...) -> Int

The number of core chemical. For instance, 2 for "[2M+H]+". 

# Generic Methods
* `AbstractAdductIon`: property search

# Specific Methods 
* Entity Level

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:ncore`
    * default: `1`
"""
ncore(adduct_ion::AbstractAdductIon; kwargs...) = getchemicalproperty(adduct_ion, :ncore, 1)

"""
    charge(chemical::AbstractChemicalsSchema; kwargs...) -> Int

The charge state of `chemical`; positive for cation and negative for anion. For instance, -2 for dianion, +3 for trication. 

# Generic Methods
* `AbstractChemicalsSchema`: property search
* `AbstractAdductIon`: charge of core chemical and adduct
* `AbstractChemicalWrapper`: charge of field `chemical`

# Specific Methods 
* Entity Level

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:charge`
    * default: `0`
"""
charge(chemical::AbstractChemicalsSchema; kwargs...) = getchemicalproperty(chemical, :charge, 0)
charge(adduct_ion::AbstractAdductIon; kwargs...) = ncore(adduct_ion) * charge(ioncore(adduct_ion); kwargs...) - charge(ionadduct(adduct_ion))
charge(chemical::AbstractChemicalWrapper; kwargs...) = charge(chemical.chemical; kwargs...)

"""
    ncharge(chemical::AbstractChemicalsSchema; kwargs...) -> Int

The number of charges of `chemical`.

# Generic Methods
* `AbstractChemicalsSchema`: absolute value of charge

# Specific Methods 
* Entity Level
"""
ncharge(chemical::AbstractChemicalsSchema; kwargs...) = abs(charge(chemical; kwargs...))

"""
    retentiontime(chemical::AbstractChemical; kwargs...) -> AbstractFloat

The retention time of `chemical`.

# Generic Methods
* `AbstractChemical`: property search
* `AbstractChemicalWrapper`: retention time of field `chemical`

# Specific Methods 
* Entity Level
* `Isobars`: weighted mean of retention times of each chemical.

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:retentiontime`
    * default: `NaN`
"""
retentiontime(chemical::AbstractChemical; kwargs...) = getchemicalproperty(chemical, :retentiontime, NaN)
retentiontime(chemical::AbstractChemicalWrapper; kwargs...) = retentiontime(chemical.chemical; kwargs...)

"""
    chemicalentity(chemical::AbstractChemical; kwargs...) -> AbstractChemical

The single chemical entity (having a single formula) from a chemical entity (i.e. itself), or chemical species. 

Attributes with Specific Methods marked as `Entity` indicate the `chemicalentity` are applied in the function for chemical types that are not inheretly single entities. 

# Generic Methods
* `AbstractChemical`: itself
* `AbstractChemicalWrapper`: wrapped chemical entity of field `chemical`

# Specific Methods
* `Isobars`: `chemicalentity(first(chemical.chemicals))`, i.e. the most abundant entity.
* `ChemicalTransition`: the very begining precursor
* `Groupedisotopomers`: the most abundant isotopomer
"""
chemicalentity(chemical::AbstractChemical; kwargs...) = chemical
chemicalentity(chemical::T; kwargs...) where {T <: AbstractChemicalWrapper} = T.name.wrapper(chemicalentity(chemical.chemical; kwargs...))

"""
    elementalscheme(scheme::AbstractScheme; kwargs...) -> ElementalSchema

The elemental scheme of `scheme`. 

# Generic Methods
* `AbstractScheme`: property search
* `AbstractStructuralScheme`: throw error

# Specific Methods
* `ChemicalSchema`: elemental schema of all schema
* `IsotopomerizedSchema`: elemental schema of all schema

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:elementalscheme`
    * default: itself
"""
elementalscheme(scheme::AbstractScheme; kwargs...) = getchemicalproperty(scheme, :elementalscheme, scheme)
elementalscheme(scheme::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("Structural scheme requires specific mapping to elemental scheme."))

"""
    elementalscheme(scheme::AbstractScheme; kwargs...) -> AbstractScheme

The structural scheme of `scheme`. All elemental schema are regarded as structural schema as well.

# Generic Methods
* `AbstractScheme`: property search

# Specific Methods
* `ChemicalSchema`: structural schema of all schema
* `IsotopomerizedSchema`: structural schema of all schema

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:structuralscheme`
    * default: itself
"""
structuralscheme(scheme::AbstractScheme; kwargs...) = getchemicalproperty(scheme, :structuralscheme, scheme)

"""
    chemicalspecies(chemical::AbstractChemical; kwargs...) -> Vector{<: AbstractChemical}

The chemical species (represented by a vector) from a chemical entity or chemical species. 

Attributes with Specific Methods marked as `Species` indicate the species are allowed with the function. 

# Generic Methods
* `AbstractChemical`: `[chemical]`

# Specific Methods
* `Isobars`: `chemical.chemicals`
* `ChemicalTransition{<: Isobars}`: a vector of chemical transition
"""
chemicalspecies(chemical::AbstractChemical; kwargs...) = [chemical]

"""
    chemicaltransition(chemical::AbstractChemical; kwargs...) -> Vector{<: AbstractChemical}

The chemical entities that are analyzed in each stage of instrumental analysis.

Attributes with Specific Methods marked as `Transition` indicate the transition are allowed with the function. 

# Generic Methods
* `AbstractChemical`: `[chemical]`

# Specific Methods
* `ChemicalTransition`: a vector of chemical entities
* `Isobars`: the most abundant chemical transition
"""
chemicaltransition(chemical::AbstractChemical; kwargs...) = [chemical]

"""
    chemicalparent(chemical::AbstractChemicalsSchema; kwargs...) -> AbstractChemical

The parent chemical without delocalized isotopes replacement.

# Generic Methods
* `AbstractChemical`: itself
* `AbstractChemicalWrapper`: wrapped parent chemical of field `chemical`
* `AbstractScheme`: itself
* `AbstractStructuralScheme`: throw error

# Specific Methods
* Entity Level
* `ChemicalTransition`: `ChemicalTransition` of parent chemicals of each transition
* `Groupedisotopmers`: field `parent` 
* `ElementalScheme`: scheme of chemical parent of chemical entity 
* `ChemicalSchema`: schema of chemical parent(s) of chemical entity(s)
* `IsotopomerizedSchema`: field `parent` 
* `Groupedisotopmerizedschema`: field `parent` 
"""
chemicalparent(cc::AbstractChemical; kwargs...) = cc
chemicalparent(chemical::T; kwargs...) where {T <: AbstractChemicalWrapper} = T.name.wrapper(chemicalparent(chemical.chemical; kwargs...))
chemicalparent(sch::AbstractScheme; kwargs...) = sch
chemicalparent(sch::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`chemicalparent` cannot be defined for structural scheme."))

"""
    isotopomersisotopes(chemical::AbstractChemicalsSchema; kwargs...) -> Vector{Pair{String, Int}}

The delocalized isotopes replacement of isotopomers.

# Generic Methods
* `AbstractChemical`: property search
* `AbstractChemicalWrapper`: delocalized isotopes replacement of field `chemical`
* `AbstractScheme`: `Pair{String, Int}[]`
* `AbstractStructuralScheme`: throw error
* `AbstractCompleteScheme`: delocalized isotopes replacement of elemental scheme

# Specific Methods 
* Entity Level
* `Isotopomers`: field `isotopes`
* `Groupedisotopomers`: isotopes replacement of the most abundunt isotopomers
* `ElementalScheme`: delocalized isotopes replacement of chemical entity loss
* `ChemicalSchema`: all delocalized isotopes replacement of chemical entity loss
* `IsotopomerizedSchema`: field `isotopes`
* `Groupedisotopomerizedschema`: isotopes replacement of the most abundunt isotopomers

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:isotopomersisotopes`
    * default: `Pair{String, Int}[]`

# Keyword arguments
* `loss` determines whether the chemical is part of chemical loss, and sign flips are propagated into elements.
"""
isotopomersisotopes(cc::AbstractChemical; kwargs...) = getchemicalproperty(cc, :isotopomersisotopes, Pair{String, Int}[])
isotopomersisotopes(chemical::AbstractChemicalWrapper; kwargs...) = isotopomersisotopes(chemical.chemical; kwargs...)
isotopomersisotopes(sch::AbstractScheme; kwargs...) = Pair{String, Int}[]
isotopomersisotopes(sch::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`isotopomersisotopes` cannot be defined for structural scheme."))
isotopomersisotopes(sch::AbstractCompleteScheme; kwargs...) = isotopomersisotopes(elementalscheme(sch); kwargs...) 

"""
    isotopomerstate(cc::AbstractChemicalsSchema; isotope = "[13C]", kwargs...) -> Int 

The isotopomers state, i.e. equivalent number of `isotope`. 

# Generic Methods
* `AbstractChemicalsSchema`: isotopomers state calculated from `isotopomersisotopes`

# Specific Methods
* Entity Level

# Keyword arguments
* `ischemical` determines whether the chemical is a chemical or a scheme. 
* `loss` determines whether the chemical is part of chemical loss, and sign flips are propagated into `isotopomersisotopes`.
"""
isotopomerstate(cc::AbstractChemical; isotope_unit = nothing, isotope = "[13C]", loss = false, ischemical = true, kwargs...) = _isotopomerstate(isotopomersisotopes(cc), isnothing(isotope_unit) ? elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]] : isotope_unit; loss, ischemical, kwargs...)
isotopomerstate(cc::AbstractScheme; isotope_unit = nothing, isotope = "[13C]", loss = true, ischemical = false, kwargs...) = _isotopomerstate(isotopomersisotopes(cc), isnothing(isotope_unit) ? elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]] : isotope_unit; loss, ischemical, kwargs...)

# truly formed chemical
"""
    analyzedchemical(chemical::AbstractChemicalsSchema; kwargs...) -> AbstractChemical

The single chemical entity that is directly analyzed in the very begining of instrumental analysis.

# Generic Methods
* `AbstractChemical`: itself
* `AbstractScheme`: throw error

# Specific Methods
* Species Level
* `ChemicalTransition`: the very begining precursor, but cannot be chemical gain or loss.
* `Isobars`: the most abundant analyzed chemical
"""
analyzedchemical(cc::AbstractChemical; kwargs...) = cc
analyzedchemical(cc::AbstractScheme; kwargs...) = throw(ArgumentError("Scheme cannot be analyzed directly."))

"""
    detectedchemical(chemical::AbstractChemicalsSchema; kwargs...) -> AbstractChemical

The single chemical entity that is directly detected in the very ending of instrumental analysis.

# Generic Methods
* `AbstractChemical`: itself
* `AbstractScheme`: precursor with scheme. It requires a keyword argument `precursor`.

# Specific Methods
* Species Level
* `ChemicalTransition`: the very ending product, but takes chemical gain and loss in consideration.
* `Isobars`: the most abundant detected chemical
"""
detectedchemical(cc::AbstractChemical; kwargs...) = cc
function detectedchemical(sch::AbstractScheme; precursor = nothing, kwargs...) 
    isnothing(precursor) && throw(ArgumentError("Scheme cannot be directly detected without precursor."))
    precursor = detectedchemical(precursor; kwargs...)
    detectedchemical(precursor, sch)
end
"""
    detectedisotopes(chemical::AbstractChemicalsSchema; kwargs...) -> Vector{Pair{String, Int}}

The delocalized isotopes replacement of detected chemical. See `detectedchemical` for details.
"""
detectedisotopes(cc::AbstractChemicalsSchema; kwargs...) = isotopomersisotopes(detectedchemical(cc); kwargs...)

"""
    detectedcharge(chemical::AbstractChemicalsSchema; kwargs...) -> Int

The charge state of detected chemical. See `detectedchemical` for details.
"""
detectedcharge(cc::AbstractChemicalsSchema; kwargs...) = charge(detectedchemical(cc); kwargs...)

"""
    detectedelements(chemical::AbstractChemicalsSchema); kwargs... -> Vector{Pair{String, Int}}

The elements of detected chemical. See `detectedchemical` for details.
"""
detectedelements(cc::AbstractChemicalsSchema; kwargs...) = chemicalelements(detectedchemical(cc); kwargs...)

# MS representation of chemical
"""
    inputchemical(chemical::AbstractChemicalsSchema; kwargs...) -> AbstractChemical

The single chemical entity that is the input of the very begining of instrumental analysis. 

It is equivalent to `analyzedchemical`, except schema are accepeted. See `analyzedchemical` for details.
"""
inputchemical(cc::AbstractChemicalsSchema; kwargs...) = cc

"""
    outputchemical(chemical::AbstractChemicalsSchema; kwargs...) -> AbstractChemicalsSchema

The single chemical entity that is the output of the very ending of instrumental analysis.

It is equivalent to `detectedchemical`, except schema are kept unchanged. See `detectedchemical` for details.
"""
outputchemical(cc::AbstractChemicalsSchema; kwargs...) = cc

"""
    seriesanalyzedchemical(chemical::AbstractChemicalsSchema; kwargs...) -> Vector{<: AbstractChemical}

The chemical entities that are directly analyzed in each stage of instrumental analysis. 

It is equivalent to `chemicaltransition`, except schema are transformed to detected chemicals.

# Generic Methods
* `AbstractChemicalsSchema`: `[chemical]`

# Specific Methods
* Transition Level
* `ChemicalTransition`: a vector of analyzed chemicals. See `analyzedchemical` for details.
"""
seriesanalyzedchemical(cc::AbstractChemicalsSchema; kwargs...) = [cc]

"""
    seriesanalyzedisotopes(chemical::AbstractChemicalsSchema; kwargs...) -> Vector{Vector{Pair{String, Int}}}

The delocalized isotopes replacement of serially analyzed chemical. See `seriesanalyzedchemical` for details.
"""
seriesanalyzedisotopes(cc::AbstractChemicalsSchema; kwargs...) = [isotopomersisotopes(c; kwargs...) for c in seriesanalyzedchemical(cc)]

"""
    seriesanalyzedcharge(chemical::AbstractChemicalsSchema; kwargs...) -> Vector{Int}

The charge states of serially analyzed chemical. See `seriesanalyzedchemical` for details.
"""
seriesanalyzedcharge(cc::AbstractChemicalsSchema; kwargs...) = [charge(c; kwargs...) for c in seriesanalyzedchemical(cc)]

"""
    seriesanalyzedelements(chemical::AbstractChemicalsSchema; kwargs...) -> Vector{Vector{Pair{String, Int}}}

The elements of serially analyzed chemical. See `seriesanalyzedchemical` for details.
"""
seriesanalyzedelements(cc::AbstractChemicalsSchema; kwargs...) = [chemicalelements(c; kwargs...) for c in seriesanalyzedchemical(cc)]

"""
    msstage(chemical::AbstractChemical; kwargs...) -> Int

Number of stages of MS the chemical has been through.

# Generic Methods
* `AbstractChemical`: 1
"""
msstage(cc::AbstractChemical; kwargs...) = 1

"""
    mmi(chemical::AbstractChemical) -> AbstractFloat

The monoisotopic mass of `chemical`.

# Generic Methods
* `AbstractChemical`: calculated from `chemicalelements` and `charge`

# Specific Methods
* Entity Level
* `Isobars`: weighted mean of monoisotopic masses
* `Groupedisotopomers`: weighted mean of monoisotopic masses
* `ChemicalTransition`: monoisotopic mass of `analyzedchemical`
"""
mmi(cc::AbstractChemical) = mmi(chemicalelements(cc), charge(cc))

"""
    molarmass(chemical::AbstractChemical) -> AbstractFloat

The molar mass of `chemical`.

# Generic Methods
* `AbstractChemical`: calculated from `chemicalelements` and `charge`

# Specific Methods
* Entity Level
* `Isobars`: weighted mean of molar masses
* `Groupedisotopomers`: weighted mean of molar masses
* `ChemicalTransition`: molar mass of `analyzedchemical`
"""
molarmass(cc::AbstractChemical) = molarmass(chemicalelements(cc), charge(cc))

"""
    mz(chemical::AbstractChemical[, adduct]) -> AbstractFloat

The mass to charge ratio (m/z) of charged chemical or chemical with adduct. It is equivalent to `mmi(charged_chemical) / ncharge(charged_chemical)`.

# Generic Methods
* `AbstractChemical`: calculated from `mmi` and `charge`
* `AbstractAdductIon`: calculated from `mmi` and `charge` of both core chemical and adduct

# Specific Methods
* Entity Level
* `Isobars`: weighted mean of m/z
* `Groupedisotopomers`: weighted mean of m/z
* `ChemicalTransition`: m/z of `analyzedchemical`
"""
mz(charged_cc::AbstractChemical) = charge(charged_cc) == 0 ? NaN : mmi(charged_cc) / ncharge(charged_cc)
mz(cc::AbstractChemical, adduct) = mz(ionize(cc; parse_adduct(adduct)...))
function mz(adduct_ion::AbstractAdductIon) 
    (ncore(adduct_ion) * mmi(ioncore(adduct_ion)) + mmi(chemicalelements(ionadduct(adduct_ion); loss = false)) - charge(adduct_ion) * ME) / ncharge(adduct_ion)
end

# list all propertys
# add new propertys
# function chemicalpropertynames() end