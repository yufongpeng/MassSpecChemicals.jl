# MassSpecChemicals

[![Build Status](https://github.com/yufongpeng/MassSpecChemicals.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/yufongpeng/MassSpecChemicals.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/yufongpeng/MassSpecChemicals.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/yufongpeng/MassSpecChemicals.jl)

A package for representing molecules or ions formed in mass spectrometers (MS). 

All chemicals are instances of abstract type `AbstractChemical`.
Charged chemicals with specific adduct or molecule loss that formed in MS (adduct ion) are instances of abstract type `AbstractAdductIon`.

# Built-in chemical types
1. `Chemical`: unstructured chemicals, storing name, formula, and other attributes; `Chemical(name::String, formula::String, attr::Vector{Pair{Symbol, Any}})`.
2. `AdductIon`: charged chemicals with specific adduct or molecule loss; `AdductIon(core::AbstractChemical, adduct::AbstractAdduct)`.
3. `Isobars`: multiple chemicals with similar m/z; `Isobars(chemical::Vector{<: AbstractChemical}, abundance::Vector{Float64})`.
4. `Isotopomers`: multiple chemicals differed from isotopic replacement location; `Isotopomers(parent::AbstractChemical, isotopes::Vector{Pair{String, Int}})`.

# Formula
## Parent Elements and Major isotopes
|Symbol|Major isotopes|Atomic number|Mass number|
|--------|-------------|-----------|------------------|
|C|[12C]|6|12|
|H|[1H]|1|1|
|O|[16H]|8|16|
|N|[14N]|7|14|
|P|[31P]|15|31|
|S|[32S]|16|32|
|Li|[7Li]|3|7|
|Na|[23Na]|11|23|
|K|[39K]|19|39|
|F|[19F]|9|19|
|Cl|[35Cl]|17|35|
|Ag|[108Ag]|47|108|

## Minor Isotopes
|Minor isotopes|Atomic number|Mass number|Alternative symbol|
|--------|-------------|-----------|------------------|
|[13C]|6|13||
|D|1|2|[2H]|
|[17O]|8|17||
|[18O]|8|18||
|[15N]|7|15||
|[33S]|16|33||
|[34S]|16|34||
|[36S]|16|36||
|[6Li]|3|6||
|[40K]|19|40||
|[41K]|19|41||
|[37Cl]|17|37||
|[109Ag]|47|109||


By default, parent elements are considered as major isotopes possibly replaced by minor isotopes. For instance,
* CO2 has a carbon-12 and two oxygen-16, but any minor isotopes replacements are possible.
* [13C][16O]O has a carbon-13, an oxygen-16, and an oxygen-16 possibly replaced by other minor isotopes.

One exception is that in `parent` chemical of `Isotopomers`, parent elements are major isotopes, and the number of replacement is restricted by field `isotopes`. 

# Adduct
All adducts are instances of abstract type `AbstractAdduct` and `AbstractPosAdduct` or `AbstractNegAdduct`.

Prefefined adducts:
|Adduct Object|Adduct expression|
|-------------|-----------------|
|`LossElectron`|[M]+|
|`Protonation`|[M+H]+|
|`ProtonationNLH2O`|[M+H-H2O]+|
|`ProtonationNL2H2O`|[M+H-2H2O]+|
|`ProtonationNL3H2O`|[M+H-3H2O]+|
|`DiProtonation`|[M+2H]2+|
|`TriProtonation`|[M+3H]3+|
|`AddNH4`|[M+NH4]+|
|`AddHNH4`|[M+H+NH4]2+|
|`Add2NH4`|[M+2NH4]2+|
|`Sodization`|[M+Na]+|
|`SodizationProtonation`|[M+Na+H]2+|
|`DiSodization`|[M+2Na]2+|
|`AddElectron`|[M]-|
|`Deprotonation`|[M-H]-|
|`DeprotonationNLH2O`|[M-H-H2O]-|
|`DiDeprotonation`|[M-2H]2-|
|`TriDeprotonation`|[M-3H]3-|
|`AddOAc`|[M+OAc]-|
|`AddFo`|[M+OFo]-|
|`Fluoridation`|[M+F]-|
|`Chloridation`|[M+Cl]-|
|`Demethylation`|[M-Me]-|

User can custumize adduct by `PosAdduct` and `NegAdduct` or define subtypes of `AbstractPosAdduct` or `AbstractNegAdduct` for more general usage.

For example, [2M+H]+, and [M-2H-2H2O]2- can be created by `PosAdduct(2, "+H", 1)`, and `NegAdduct(1, "-2H-2H2O", 2)` respectively.

# API
## Attributes of `AbstractChemical`
`getchemicalattr(chemical::ChemicalType, attr::Symbol; kwargs...)` returns the attribute `attr` of `chemical` with additional settings by `kwargs`. To define new attributes or overwrite defaults, define the following method, `getchemicalattr(chemical::ChemicalType, attr::Val{attr}; kwargs...)`.

The following getter functions are user friendly interfaces for `getchemicalattr`.
* `chemicalname`: unique chemical name (`String`).
* `chemicalformula`: chemical formula (`String`). 
* `chemicalelements`: chemical elements (`Vector{Pair{String, Int}}`). 
* `chemicalabbr`: common abbreviation (`String`); defaults to `chemicalname`.  
* `chemicalsmiles`: SMILES (`String`); defaults to `""`. 
* `charge`: charges (positive or negative); defaults to 0. 
* `abundant_chemical`: the most abundant chemical from a chemical (itself) or isobars.
* `rt`: retention time; defaults `NaN`.  

### Additional Atttributes of `AbstractAdductIon`
* `ioncore`: core of an adduct ion.
* `ionadduct`: adduct of an adduct ion.
* `kmer`: number of core chemical "M".

The following function should be extended when isotopes are involved in addut ion formation.
* `adductisotopes`: element-number pairs; the elements changed when the core chemical has isotopic labeling that is lost in adduct formation. For instance, [M-Me]- of D-labeled PC may turn out to be [M-CD3]- rather than [M-CH3]- if Ds are on the methyl group. In this case, `adductisotopes` of [M-Me]- of PC should be `["H" => 3, "D" => -3]`.

### Additional Atttributes of `Isobars`
* `chemicals`: composition of isabars.
* `abundance`: abundance of each isabars.

### Additional Atttributes of `Isotopomers`
* `parent`: shared chemical structure of isotopomers prior to isotopic replacement.
* `isotopes`: isotopes-number pairs of isotopic replacement.

### Additional Atttributes of `ChemicalPair`
* `precursor`: chemical prior to fragmentation.
* `product`: fragments.

### Functions for Derived Attributes
* `mmi`: monoisotopic mass.
* `molarmass`: molar mass.
* `mz`: m/z.
* `ncharge`: number of charges. 

## Functions for Attributes of `AbstractAdduct`
* `kmer`: number of core chemical "M".
* `charge`: charges.
* `ncharge`: number of charges. 
* `adductformula`: the formula for adduct. For instance, `"-H"` for [M-H]-, `"+OAc"` for [M+OAc]-.
* `adductelements`: elements changed with adduct formation. For generic adduct, it uses `adductformula` to calculate elements. If non-element strings are used in `adductformula`, defining custumized `adductelements` is required. 

## Other Functions
* `ischemicalequal`: whether two chemicals equal.
* `isadductequal`: whether two adducts equal.
* `isotopicabundance`: isotopic abundance of a single chemical.
* `isotopologues`: isotopologues of a single chemical or formula, and MS/MS transition.
* `isotopologues_table`: a table with isotopologues, and more information. 
* `isobars_rt`: isobars under specific m/z, mw, and rt tolerance.
* `isobars_rt_table`: a table with isobars, and more information. 
* `acrit`: create absolute criterion.
* `rcrit`: create relative criterion.
* `crit`: create both absolute and relative criterion.
* `@ri_str`: real number interval.
