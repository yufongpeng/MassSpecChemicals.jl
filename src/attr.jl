defaultname(::T) where {T <: AbstractChemical} = string("Chemical::", T)

"""
    getchemicalproperty(chemical::AbstractChemical, property::Symbol, default = nothing)

Get property from `chemical`. 

This function defaults to finds the property, and returns `default` If it is not available.

For type `Chemical` and properties other than `name`, `elements` and `formula`, it iterates through `chemical.property`. If no matched property name is found, it returns `default`.

For type `AbstractAdductIon`, it searches for properties of itself and then the properties of core chemical without specialized methods.
"""
getchemicalproperty(chemical::AbstractChemical, property::Symbol, default = nothing) = hasproperty(chemical, property) ? getproperty(chemical, property) : default 
getchemicalproperty(adduct_ion::AbstractAdductIon, property::Symbol, default) = hasproperty(adduct_ion, property) ? getproperty(adduct_ion, property) : getchemicalproperty(ioncore(adduct_ion), property, default)

function getchemicalproperty(chemical::Union{Chemical, FormulaChemical}, property::Symbol, default)
    hasproperty(chemical, property) && return getproperty(chemical, property)
    for (p, v) in chemical.property
        p == property && return v
    end
    return default
end

"""
    chemicalname(chemical::AbstractChemical; verbose = true, kwargs...) -> String

The name of `chemical`. It should be unique for each chemical. 

# Generic Methods
* `AbstractChemical`: property search
* `AbstractAdductIon`: name composed of name of core chemical and adduct

# Specific Methods 
* Species/Transition Level

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:name`
    * default: `defaultname(chemical)`

# Keyword arguments
* `verbose` determines whether includes all names or not for chemical species. If `verbose` is false, only the first (most abundant) chemical is included.
"""
chemicalname(chemical::AbstractChemical; kwargs...) = getchemicalproperty(chemical, :name, defaultname(chemical))
function chemicalname(adduct_ion::AbstractAdductIon; corename = nothing, kwargs...) 
    r = isnothing(corename) ? chemicalname(ioncore(adduct_ion); kwargs...) : corename
    if occursin(" ", r)
        r = string("(", r, ")")
    end
    if isnothing(corename)
        s = replace(string(ionadduct(adduct_ion)), "M" => r; count = 1)
        s = replace(s, r"\d*[+-]$" => ""; count = 1)
    else
        s = replace(string(ionadduct(adduct_ion)), r"\d*M" => r; count = 1)
        s = replace(s, r"\d*[+-]$" => ""; count = 1)
    end
    c = charge(adduct_ion)
    c == 0 ? s : string(s, abs(c) > 1 ? abs(c) : "", c > 0 ? "+" : "-") 
end

"""
    chemicalformula(chemical::AbstractChemical; kwargs...) -> String

The formula of chemical entity of `chemical`.

# Generic Methods
* `AbstractChemical`: property search
* `AbstractAdductIon`: formuala combining core chemical and adduct

# Specific Methods 
* Entity Level

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

function chemicalformula(adduct_ion::AbstractAdductIon; kwargs...)
    el = deepcopy(chemicalelements(ioncore(adduct_ion); kwargs...))
    nm = kmer(adduct_ion)
    if nm > 1
        for id in eachindex(el)
            el[id] = first(el[id]) => nm * last(el[id])
        end
    end
    chemicalformula(gain_elements!(dictionary_elements(el), adductelements(adduct_ion)))
end

"""
    chemicalelements(chemical::AbstractChemical; kwargs...) -> Vector{Pair{String, Int}}

The elements of chemical entity of `chemical`.

# Generic Methods
* `AbstractChemical`: property search
* `AbstractAdductIon`: elements combining core chemical and adduct

# Specific Methods 
* Entity Level

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

function chemicalelements(adduct_ion::AbstractAdductIon; kwargs...) 
    el = deepcopy(chemicalelements(ioncore(adduct_ion); kwargs...))
    nm = kmer(adduct_ion)
    if nm > 1
        for id in eachindex(el)
            el[id] = first(el[id]) => nm * last(el[id])
        end
    end
    vcat(el, adductelements(adduct_ion))
end

"""
    chemicalabbr(chemical::AbstractChemical; verbose = true, kwargs...) -> String

The abbreviation of `chemical`. 

# Target Chemical Level
* Species

# Generic Methods
* `AbstractChemical`: property search
* `AbstractAdductIon`: abbreviation composed of abbreviation of core chemical and adduct

# Specific Methods 
* Species/Transition Level

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:abbreviation`
    * default: `chemicalname(chemical; kwargs...)`

# Keyword arguments
* `verbose` determines whether includes all abbreviations or not for chemical species. If `verbose` is false, only the first (most abundant) chemical is included.
"""
chemicalabbr(chemical::AbstractChemical; kwargs...) = getchemicalproperty(chemical, :abbreviation, chemicalname(chemical; kwargs...))
function chemicalabbr(adduct_ion::AbstractAdductIon; kwargs...) 
    r = chemicalabbr(ioncore(adduct_ion); kwargs...)
    if occursin(" ", r)
        r = string("(", r, ")")
    end
    s = replace(string(ionadduct(adduct_ion)), "M" => r, r"\d*[+-]$" => ""; count = 2)
    c = charge(adduct_ion)
    c == 0 ? s : string(s, abs(c) > 1 ? abs(c) : "", c > 0 ? "+" : "-") 
end

"""
    chemicalsmiles(chemical::AbstractChemical; kwargs...) -> String

The SMILES of chemical entity of `chemical`.

# Generic Methods
* `AbstractChemical`: property search
* `AbstractAdductIon`: SMILES of the core chemical

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
    kmer(adduct_ion::AbstractAdductIon; kwargs...) -> Int

The number of core chemical. For instance, 2 for "[2M+H]+". 

# Generic Methods
* `AbstractAdductIon`: `kmer` of `ionadduct`

# Specific Methods 
* Entity Level
"""
kmer(adduct_ion::AbstractAdductIon; kwargs...) = kmer(ionadduct(adduct_ion))

"""
    charge(chemical::AbstractChemical; kwargs...) -> Int

The charge state of `chemical`; positive for cation and negative for anion. For instance, -2 for dianion, +3 for trication. 

# Generic Methods
* `AbstractChemical`: property search
* `AbstractAdductIon`: charge of core chemical and adduct

# Specific Methods 
* Entity Level

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

# Generic Methods
* `AbstractChemical`: absolute value of charge

# Specific Methods 
* Entity Level
"""
ncharge(chemical::AbstractChemical; kwargs...) = abs(charge(chemical; kwargs...))

"""
    retentiontime(chemical::AbstractChemical; kwargs...) -> AbstractFloat

The retention time of `chemical`.

# Generic Methods
* `AbstractChemical`: property search

# Specific Methods 
* Entity Level
* `Isobars`: weighted mean of retention times of each chemical.

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:retentiontime`
    * default: `NaN`
"""
retentiontime(chemical::AbstractChemical; kwargs...) = getchemicalproperty(chemical, :retentiontime, NaN)

"""
    chemicalentity(chemical::AbstractChemical) -> AbstractChemical

The single chemical entity (having a single formula) from a chemical entity (i.e. itself), or chemical species. 

Attributes with Specific Methods marked as `Entity` indicate the `chemicalentity` are applied in the function for chemical types that are not inheretly single entities. 

# Generic Methods
* `AbstractChemical`: itself

# Specific Methods
* `Isobars`: `chemicalentity(first(chemical.chemicals))`, i.e. the most abundant entity.
* `ChemicalLoss`: `chemical.chemical`
* `ChemicalGain`: `chemical.chemical`
* `ChemicalTransition`: the very begining precursor
* `Groupedisotopomers`: the most abundant isotopomer
"""
chemicalentity(chemical::AbstractChemical; kwargs...) = chemical

"""
    chemicalspecies(chemical::AbstractChemical) -> Vector{<: AbstractChemical}

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
    chemicaltransition(chemical::AbstractChemical) -> Vector{<: AbstractChemical}

The chemical entities that are analyzed in each stage of instrumental analysis.

Attributes with Specific Methods marked as `Transition` indicate the transition are allowed with the function. 

# Generic Methods
* `AbstractChemical`: `[chemical]`

# Specific Methods
* `ChemicalTransition`: a vector of chemical entities
* `Isobars`: the most abundant chemical transition
"""
chemicaltransition(chemical::AbstractChemical) = [chemical]

"""
    chemicalparent(chemical::AbstractChemical) -> AbstractChemical

The parent chemical without delocalized isotopes replacement.

# Generic Methods
* `AbstractChemical`: itself

# Specific Methods
* Entity Level
*`ChemicalTransition`: `ChemicalTransition` of parent chemicals of each transition
"""
chemicalparent(cc::AbstractChemical; kwargs...) = cc

"""
    isotopomersisotopes(chemical::AbstractChemical) -> Vector{Pair{String, Int}}

The delocalized isotopes replacement of isotopomers.

# Generic Methods
* `AbstractChemical`: property search

# Specific Methods 
* Entity Level
* `Isotopomers`: `chemical.isotopes`
* `Groupedisotopomers`: isotopes replacement of the most abundunt isotopomers

# Property Search Workflow
1. Function: `getchemicalproperty`
    * property: `:isotopomersisotopes`
    * default: `Pair{String, Int}[]`
"""
isotopomersisotopes(cc::AbstractChemical; kwargs...) = getchemicalproperty(cc, :isotopomersisotopes, Pair{String, Int}[])

"""
    isotopomerstate(cc::AbstractChemical; isotope = "[13C]") -> Int 

The isotopomers state, i.e. equivalent number of `isotope`. 

# Generic Methods
* `AbstractChemical`: isotopomers state calculated from `isotopomersisotopes`

# Specific Methods
* Entity Level
"""
isotopomerstate(cc::AbstractChemical; isotope = "[13C]") = _isotopomerstate(isotopomersisotopes(cc), elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]])

# truly formed chemical
"""
    analyzedchemical(chemical::AbstractChemical) -> AbstractChemical

The single chemical entity that is directly analyzed in the very begining of instrumental analysis.

# Generic Methods
* `AbstractChemical`: itself

# Specific Methods
* Species Level
* `ChemicalLoss`: throw error
* `ChemicalGain`: throw error
* `ChemicalTransition`: the very begining precursor, but cannot be chemical gain or loss.
* `Isobars`: the most abundant analyzed chemical
"""
analyzedchemical(cc::AbstractChemical; kwargs...) = cc

"""
    detectedchemical(chemical::AbstractChemical) -> AbstractChemical

The single chemical entity that is directly detected in the very ending of instrumental analysis.

# Generic Methods
* `AbstractChemical`: itself

# Specific Methods
* Species Level
* `ChemicalLoss`: precursor exclueded `chemical.chemical` part. It requires a keyword argument `precursor`.
* `ChemicalGain`: precursor inclueded `chemical.chemical` part. It requires a keyword argument `precursor`.
* `ChemicalTransition`: the very ending product, but takes chemical gain and loss in consideration.
* `Isobars`: the most abundant detected chemical
"""
detectedchemical(cc::AbstractChemical; kwargs...) = cc

"""
    detectedisotopes(chemical::AbstractChemical) -> Vector{Pair{String, Int}}

The delocalized isotopes replacement of detected chemical. See `detectedchemical` for details.
"""
detectedisotopes(cc::AbstractChemical; kwargs...) = isotopomersisotopes(detectedchemical(cc); kwargs...)

"""
    detectedcharge(chemical::AbstractChemical) -> Int

The charge state of detected chemical. See `detectedchemical` for details.
"""
detectedcharge(cc::AbstractChemical; kwargs...) = charge(detectedchemical(cc); kwargs...)

"""
    detectedelements(chemical::AbstractChemical) -> Vector{Pair{String, Int}}

The elements of detected chemical. See `detectedchemical` for details.
"""
detectedelements(cc::AbstractChemical; kwargs...) = chemicalelements(detectedchemical(cc); kwargs...)

# MS representation of chemical
"""
    inputchemical(chemical::AbstractChemical) -> AbstractChemical

The single chemical entity that is the input of the very begining of instrumental analysis. 

It is equivalent to `analyzedchemical`, except `ChemicalLoss` and `ChemicalGain` are accepeted. See `analyzedchemical` for details.
"""
inputchemical(cc::AbstractChemical; kwargs...) = cc

"""
    outputchemical(chemical::AbstractChemical) -> AbstractChemical

The single chemical entity that is the output of the very ending of instrumental analysis.

It is equivalent to `detectedchemical`, except `ChemicalLoss` and `ChemicalGain` are kept unchanged. See `detectedchemical` for details.
"""
outputchemical(cc::AbstractChemical; kwargs...) = cc

"""
    seriesanalyzedchemical(chemical::AbstractChemical) -> Vector{<: AbstractChemical}

The chemical entities that are directly analyzed in each stage of instrumental analysis. 

It is equivalent to `chemicaltransition`, except that `ChemicalLoss` and `ChemicalGain` are transformed to detected chemicals.

# Generic Methods
* `AbstractChemical`: `[chemical]`

# Specific Methods
* Transition Level
* `ChemicalTransition`: a vector of analyzed chemicals. See `analyzedchemical` for details.
"""
seriesanalyzedchemical(cc::AbstractChemical; kwargs...) = [cc]

"""
    seriesanalyzedisotopes(chemical::AbstractChemical) -> Vector{Vector{Pair{String, Int}}}

The delocalized isotopes replacement of serially analyzed chemical. See `seriesanalyzedchemical` for details.
"""
seriesanalyzedisotopes(cc::AbstractChemical; kwargs...) = [isotopomersisotopes(c; kwargs...) for c in seriesanalyzedchemical(cc)]

"""
    seriesanalyzedcharge(chemical::AbstractChemical) -> Vector{Int}

The charge states of serially analyzed chemical. See `seriesanalyzedchemical` for details.
"""
seriesanalyzedcharge(cc::AbstractChemical; kwargs...) = [charge(c; kwargs...) for c in seriesanalyzedchemical(cc)]

"""
    seriesanalyzedelements(chemical::AbstractChemical) -> Vector{Vector{Pair{String, Int}}}

The elements of serially analyzed chemical. See `seriesanalyzedchemical` for details.
"""
seriesanalyzedelements(cc::AbstractChemical; kwargs...) = [chemicalelements(c; kwargs...) for c in seriesanalyzedchemical(cc)]

"""
    msstage(chemical::AbstractChemical) -> Int

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
mz(cc::AbstractChemical, adduct) = mz(AdductIon(cc, parse_adduct(adduct)))
function mz(adduct_ion::AbstractAdductIon) 
    adduct = ionadduct(adduct_ion)
    (kmer(adduct) * mmi(ioncore(adduct_ion)) + mmi(adductelements(adduct_ion)) - charge(adduct_ion) * ustrip(ME)) / ncharge(adduct_ion)
end
mz(adduct_ion::AbstractAdductIon, adduct) = mz(AdductIon(ioncore(adduct_ion), parse_adduct(adduct)))

"""
    adductelements(adduct_ion::AbstractAdductIon) -> Vector{Pair{String, Int}}

The elements changed with adduct of `adduct_ion`. It contains the elements of adduct itself and isotopic labeling related to the adduct for the `adduct_ion` (`adductisotopes(adduct_ion)`). 

# Generic Methods
* `AbstractAdductIon`: calculated from `adductelements` of core chemical and `adductisotopes` of the adduct ion
"""
adductelements(adduct_ion::AbstractAdductIon) = vcat(adductelements(ionadduct(adduct_ion)), adductisotopes(adduct_ion))

"""
    adductisotopes(adduct_ion::AbstractAdductIon) -> Vector{Pair{String, Int}}

The elements changed when the core chemical has isotopic labeling that is lost in adduct formation. The returned vector is element-number pairs.

For instance, [M-Me]- (`Demethylation`) of Deuterium-labeled phosphatidylcholine (PC) may turn out to be [M-CD3]- rather than [M-CH3]- if Deuteriums are labeled on the methyl group of choline (`DLMC_PC`). In this case, `adductisotopes(::AbstractAdductIon{Demethylation, DLMC_PC})` should return `["H" => 3, "D" => -3]`.

# Generic Methods
* `AbstractAdductIon`: `Pair{String, Int}[]`

# Specific Methods
* `AbstractAdductIon{Chemical}`: user can define an additional property `:adductisotopes` for core chemical. The property should be ionadduct-(elements-number pairs) pairs. 
This function finds this property, and extracts the value of key `ionadduct(adduct_ion)`. If the property or the key does not exist, empty vector is returned. 
"""
adductisotopes(adduct_ion::AbstractAdductIon) = Pair{String, Int}[] # ex 1D: ["H" => 1, "D" => -1]
function adductisotopes(adduct_ion::AbstractAdductIon{Chemical})
    a = ionadduct(adduct_ion)
    v = getchemicalproperty(ioncore(adduct_ion), :adductisotopes, Pair{String, Int}[])
    isempty(v) && return v
    i = findfirst(x -> first(x) == a, v)
    isnothing(i) ? Pair{String, Int}[] : convert(Vector{Pair{String, Int}}, last(v[i])) # Force convert to pair as property does not restrict input type
end

# list all propertys
# add new propertys
# function chemicalpropertynames() end