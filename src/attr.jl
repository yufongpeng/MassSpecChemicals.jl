"""
    getchemicalattr(chemical::AbstractChemicalScheme) -> Union{ChemicalAttribute, Nothing}

Get where the `chemical` stores additional attributes. The name of storage depends on `getchemicalattrname`.
"""
function getchemicalattr(chemical::T) where {T <: AbstractChemicalScheme} 
    attr = getchemicalattrname(T)
    isnothing(attr) ? nothing : getproperty(chemical, attr)
end

"""
    getchemicalattrname(chemical::AbstractChemicalScheme) -> Union{Symbol, Nothing}
    getchemicalattrname(::Type) -> Union{Symbol, Nothing}

Get the name of storage where `chemical` stores additional attributes.
`AttributedChemical` has default `:attribute`, while others have no default name.
"""
getchemicalattrname(::T) where {T <: AbstractChemicalScheme} = getchemicalattrname(T)
getchemicalattrname(::Type{<: AttributedChemical}) = :attribute
getchemicalattrname(::Type{<: AbstractChemicalScheme}) = nothing

"""
    pushchemicalattr!(chemical::AbstractChemicalScheme, attr::Symbol, value)
    pushchemicalattr!(chemical::AbstractChemicalScheme, attr::AbstractString, value)

Push `value` into `chemical` as a new attribute `attr`. If `attr` already exists, the `value` replaces the old one.
"""
pushchemicalattr!(chemical::AbstractChemicalScheme, attr::AbstractString, value) = pushchemicalattr!(chemical, Symbol(attr), value)
function pushchemicalattr!(chemical::T, attr::Symbol, value) where {T <: AbstractChemicalScheme}
    hasproperty(T, attr) && return setproperty!(chemical, attr, value)
    attribute = getchemicalattr(chemical)
    isnothing(attribute) && return pushchemicalattr_new!(chemical, attr, value)
    if isnothing(attribute.attribute)
        attribute.attribute = Pair{Symbol, Any}[attr => value]
    else
        id = findfirst(==(attr) ∘ first, attribute.attribute)
        if isnothing(id)
            push!(attribute.attribute, attr => value)
        else
            attribute.attribute[id] = attr => value
        end
    end
end

function pushchemicalattr_new!(chemical::T, attr::Symbol, value) where {T <: AbstractChemicalScheme}
    setfield!(chemical, getchemicalattrname(T), ChemicalAttribute(Pair{Symbol, Any}[attr => value]))
end

"""
    getchemicalattr(chemical::AbstractChemicalScheme, attr::Symbol, default::Nothing, type = Any) -> Any
    getchemicalattr(chemical::AbstractChemicalScheme, attr::Symbol, default::T, type = T) -> T

Fetch `attr` from `chemical`.

This function looks up `attr` and returns `default` when it is unavailable. The returned value is asserted to be of type `type`.

It checks the chemical's properties before searching `chemicalattr(chemical)`;
for [`AbstractAdductIon`](@ref), it checks the adduct itself first and then falls back to the core chemical when no specialized method is available;
If no matching name is found, it returns `default`.
"""
getchemicalattr(chemical::AbstractChemicalScheme, attr::Symbol, default::T) where T = getchemicalattr(chemical, attr, default, T)
getchemicalattr(chemical::AbstractChemicalScheme, attr::Symbol, default::Nothing) = getchemicalattr(chemical, attr, default, Any)
# getchemicalattr(chemical::AbstractChemicalScheme, attr::Symbol, default, ::Type{T}) where T = hasproperty(chemical, attr) ? getproperty(chemical, attr)::T : default::T 
getchemicalattr(adduct_ion::AbstractAdductIon, attr::Symbol, default, type::Type{T})  where T = hasproperty(adduct_ion, attr) ? getproperty(adduct_ion, attr)::T : getchemicalattr(ioncore(adduct_ion), attr, default, type)
function getchemicalattr(chemical::AbstractChemicalScheme, attr::Symbol, default, ::Type{T}) where T
    hasproperty(chemical, attr) && return getproperty(chemical, attr)::T
    attribute = getchemicalattr(chemical)
    (isnothing(attribute) || isnothing(attribute.attribute)) && return default::T
    for (p, v) in attribute.attribute
        p == attr && return v::T
    end
    default::T
end

"""
    chemicalname(chemical::AbstractChemicalScheme; verbose = true, n = 1, loss = false, bracket = true, kwargs...) -> String

The name of `chemical`. This should be unique for each chemical object.

# Generic Methods
* `AbstractChemicalScheme`: attribute seach.
* `AbstractAdductIon`: name composed of name of core chemical and adduct.

# Specific Methods
* Species/Transition Level.

# Attribute Seach Workflow
1. Function: [`getchemicalattr`](@ref).
    * attribute: `:name`.
    * default: [`defaultname`](@ref).

# Keyword arguments
* `verbose::Bool` selects whether to include all names for a chemical species. If `verbose` is false, only the first, most abundant chemical is included.
* `n::Int` selects the chemical count.
* `loss::Bool` determines whether the chemical is part of a loss process, with sign changes propagated to the chemical.
* `bracket::Bool` determines whether the output string is wrapped in "[...]" with the charge included.
"""
function chemicalname(chemical::AbstractChemicalScheme; n = 1, kwargs...) 
    nm = getchemicalattr(chemical, :name, defaultname(chemical))
    n > 1 ? string(n, nm) : nm
end
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
    chemicalformula(chemical::AbstractChemical; delim = "", unique = true, ischemical = true, loss = false, kwargs...) -> String
    chemicalformula(chemical::AbstractScheme; delim = "", unique = true, ischemical = false, loss = false, kwargs...) -> String

The formula of the chemical entity represented by `chemical`.

# Generic Methods
* `AbstractChemical`: attribute seach.
* `AbstractAdductIon`: formula combining the core chemical and the adduct.
* `AbstractScheme`: attribute seach; started with `"+"` for chemical gain and `"-"` for chemical loss.
* `AbstractStructuralScheme`: throw error.

# Specific Methods
* Entity Level.

# Attribute Seach Workflow
1. Function: [`getchemicalattr`](@ref) -> `reverse_formula`.
    * attribute: `:formula`.
    * default: `""`.
2. Function: [`chemicalelements`](@ref) -> [`chemicalformula`](@ref).

# Keyword arguments
* `delim::Union{String, Char}` assigns the delimiter between each element.
* `unique::Bool` determines whether elements are combined into a unique representation when constructing the formula from [`chemicalelements`](@ref).
It defaults to false for [`Chemical`](@ref) and [`FormulaChemical`](@ref).
* `ischemical::Bool` determines whether the input represents a chemical or a scheme.
* `loss::Bool` determines whether the species is part of a loss process, with signs factored out from the element counts.
"""
function chemicalformula(chemical::AbstractChemical; shallow = false, loss = false, ischemical = true, kwargs...) 
    result = getchemicalattr(chemical, :formula, "")
    if isempty(result) && !shallow
        chemicalformula(chemicalelements(chemical; loss = false, shallow = true); loss, ischemical, kwargs...)
    else 
        reverse_formula(result, ischemical, loss)
    end
end

function chemicalformula(chemical::AbstractScheme; shallow = false, loss = false, ischemical = false, kwargs...) 
    result = getchemicalattr(chemical, :formula, "")
    if isempty(result) && !shallow
        chemicalformula(chemicalelements(chemical; loss = false, shallow = true); loss, ischemical, kwargs...)
    else 
        reverse_formula(result, ischemical, loss)
    end
end

function chemicalformula(adduct_ion::AbstractAdductIon; loss = false, kwargs...)
    el = copy(chemicalelements(ioncore(adduct_ion); loss, kwargs...))
    nm = ncore(adduct_ion)
    if nm > 1
        for id in eachindex(el)
            el[id] = first(el[id]) => nm * last(el[id])
        end
    end
    chemicalformula(gain_elements!(dictionary_elements(Dictionary, el), chemicalelements(ionadduct(adduct_ion); loss, kwargs...)); loss = !loss, kwargs...)
end

chemicalformula(sch::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`chemicalformula` cannot be defined for structural scheme."))

"""
    chemicalelements(chemical::AbstractChemicalScheme; loss = false, kwargs...) -> Vector{Pair{String, Int}}

The elemental composition of `chemical`.

# Generic Methods
* `AbstractChemical`: attribute seach.
* `AbstractAdductIon`: combined elements of the core chemical and the adduct.
* `AbstractScheme`: attribute seach; element changes are applied.
* `AbstractStructuralScheme`: throw error.

# Specific Methods
* Entity Level.

# Attribute Seach Workflow
1. Function: [`getchemicalattr`](@ref) -> `reverse_elements`.
    * attribute: `:elements`.
    * default: `Pair{String, Int}[]`.
2. Function: [`chemicalformula`](@ref) -> [`chemicalelements`](@ref).

# Keyword arguments
* `loss::Bool` determines whether the chemical participates in a loss process, with sign flips propagated into the element counts.
"""
function chemicalelements(chemical::AbstractChemical; shallow = false, loss = false, kwargs...) 
    result = getchemicalattr(chemical, :elements, Pair{String, Int}[])
    if isempty(result) && !shallow
        chemicalelements(chemicalformula(chemical; loss = false, shallow = true); loss, kwargs...)
    else 
        reverse_elements(result, loss)
    end
end

function chemicalelements(chemical::AbstractScheme; shallow = false, loss = false, kwargs...) 
    result = getchemicalattr(chemical, :elements, Pair{String, Int}[])
    if isempty(result) && !shallow
        chemicalelements(chemicalformula(chemical; loss = false, ischemical = false, shallow = true); loss, kwargs...)
    else 
        reverse_elements(result, loss) 
    end
end

function chemicalelements(adduct_ion::AbstractAdductIon; loss = false, kwargs...) 
    el = copy(chemicalelements(ioncore(adduct_ion); loss, kwargs...))
    nm = ncore(adduct_ion)
    if nm > 1
        for id in eachindex(el)
            el[id] = first(el[id]) => nm * last(el[id])
        end
    end
    vcat(el, chemicalelements(ionadduct(adduct_ion); loss, kwargs...))
end

chemicalelements(sch::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`chemicalelements` cannot be defined for structural scheme."))

"""
    chemicalabbr(chemical::AbstractChemicalScheme; verbose = true, n = 1, loss = false, bracket = true, kwargs...) -> String

The abbreviation of `chemical`. 

# Generic Methods
* `AbstractChemicalScheme`: attribute seach.
* `AbstractAdductIon`: abbreviation composed of abbreviation of core chemical and adduct.

# Specific Methods
* Species/Transition Level.

# Attribute Seach Workflow
1. Function: [`getchemicalattr`](@ref).
    * attribute: `:abbreviation`.
    * default: [`chemicalname`](@ref).

# Keyword arguments
* `verbose::Bool` determines whether includes all abbreviations or not for chemical species. If `verbose` is false, only the first (most abundant) chemical is included.
* `n::Int` selects the chemical count.
* `loss::Bool` determines whether the chemical is part of a loss process, with sign changes propagated to the chemical.
* `bracket::Bool` determines whether the output string is wrapped in "[...]" with the charge included.
"""
function chemicalabbr(chemical::AbstractChemicalScheme; n = 1, kwargs...) 
    abbr = getchemicalattr(chemical, :abbreviation, "")
    isempty(abbr) && return chemicalname(chemical; n, kwargs...)
    string(n > 1 ? n : "", abbr)
end
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

The SMILES representation of the chemical entity represented by `chemical`.

# Generic Methods
* `AbstractChemical`: attribute seach.
* `AbstractAdductIon`: SMILES of the core chemical.

# Specific Methods
* Entity Level.
* `Isotopomers`: SMILES of the parent chemical.
* `Groupedisotopomers`: SMILES of the most abundant parent chemical.

# Attribute Seach Workflow
1. Function: [`getchemicalattr`](@ref).
    * attribute: `:SMILES`.
    * default: `""`.
"""
chemicalsmiles(chemical::AbstractChemical; kwargs...) = getchemicalattr(chemical, :SMILES, "")
chemicalsmiles(adduct_ion::AbstractAdductIon; kwargs...) = chemicalsmiles(ioncore(adduct_ion); kwargs...) 

"""
    ioncore(adduct_ion::AbstractAdductIon{S, T}; kwargs...) -> S

The core chemical of `adduct_ion`. 

# Generic Methods
* `AbstractAdductIon`: attribute seach.

# Specific Methods
* Entity Level.

# Attribute Seach Workflow
1. Function: [`getchemicalattr`](@ref).
    * attribute: `:core`.
    * default: `TypeError`.
"""
ioncore(adduct_ion::AbstractAdductIon{S, T}; kwargs...) where {S, T} = getchemicalattr(adduct_ion, :core, nothing)::S

"""
    ionadduct(adduct_ion::AbstractAdductIon{S, T}; kwargs...) -> T

The adduct of `adduct_ion`. 

# Generic Methods
* `AbstractAdductIon`: attribute seach.

# Specific Methods
* Entity Level.

# Attribute Seach Workflow
1. Function: [`getchemicalattr`](@ref).
    * attribute: `:adduct`.
    * default: `TypeError`.
"""
ionadduct(adduct_ion::AbstractAdductIon{S, T}; kwargs...) where {S, T} = getchemicalattr(adduct_ion, :adduct, nothing)::T

"""
    ncore(adduct_ion::AbstractAdductIon; kwargs...) -> Int

The number of core chemicals. For example, `2` for "[2M+H]+".

# Generic Methods
* `AbstractAdductIon`: attribute seach.

# Specific Methods
* Entity Level.

# Attribute Seach Workflow
1. Function: [`getchemicalattr`](@ref).
    * attribute: `:ncore`.
    * default: `1`.
"""
ncore(adduct_ion::AbstractAdductIon; kwargs...) = getchemicalattr(adduct_ion, :ncore, 1)

"""
    charge(chemical::AbstractChemicalScheme; loss = false, kwargs...) -> Int

The charge state of `chemical`; positive for a cation and negative for an anion. For example, `-2` for a dianion and `+3` for a trication.

# Generic Methods
* `AbstractChemicalScheme`: attribute seach.
* `AbstractAdductIon`: charge of core chemical and adduct.
* `AbstractScheme`: attribute seach; change of charges.
* `AbstractCompleteScheme`: charge of the elemental scheme.
* `AbstractStructuralScheme`: throw error.

# Specific Methods
* Entity Level.

# Attribute Seach Workflow
1. Function: [`getchemicalattr`](@ref) -> sign flip if `loss`.
    * attribute: `:charge`.
    * default: `0`.

# Keyword arguments
* `loss::Bool` determines whether the chemical participates in a loss process, with sign flips propagated through the charge.
"""
charge(chemical::AbstractChemicalScheme; loss = false, kwargs...) = loss ? -getchemicalattr(chemical, :charge, 0) : getchemicalattr(chemical, :charge, 0)
charge(adduct_ion::AbstractAdductIon; loss = false, kwargs...) = ncore(adduct_ion) * charge(ioncore(adduct_ion); loss, kwargs...) + charge(ionadduct(adduct_ion); loss)
charge(sch::AbstractCompleteScheme; loss = false, kwargs...) = charge(elementalscheme(sch); loss, kwargs...)
charge(sch::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`charge` cannot be defined for structural scheme."))

"""
    ncharge(chemical::AbstractChemicalScheme; kwargs...) -> Int

The number of charges of `chemical`.

# Generic Methods
* `AbstractChemicalScheme`: absolute value of charge.

# Specific Methods
* Entity Level.
"""
ncharge(chemical::AbstractChemicalScheme; kwargs...) = abs(charge(chemical; kwargs...))

"""
    retentiontime(chemical::AbstractChemical; kwargs...) -> AbstractFloat

The retention time of `chemical`.

# Generic Methods
* `AbstractChemical`: attribute seach.

# Specific Methods
* Entity Level.
* `Isobars`: weighted mean of retention times of each chemical.

# Attribute Seach Workflow
1. Function: [`getchemicalattr`](@ref).
    * attribute: `:retentiontime`.
    * default: `NaN`.
"""
retentiontime(chemical::AbstractChemical; kwargs...) = getchemicalattr(chemical, :retentiontime, NaN)

"""
    chemicalentity(chemical::AbstractChemical; kwargs...) -> AbstractChemical

The single chemical entity (with a single formula) derived from a chemical entity or chemical species.

A `Specific Methods` entry marked as `Entity` indicates that `chemicalentity` is applied for a chemical type that is not originally a single entity.

# Generic Methods
* `AbstractChemical`: itself.

# Specific Methods
* `Isobars`: `chemicalentity(first(chemical.chemicals))`, i.e. the most abundant entity.
* `ChemicalTransition`: the very beginning precursor.
* `Groupedisotopomers`: the most abundant isotopomer.
"""
chemicalentity(chemical::AbstractChemical; kwargs...) = chemical

"""
    elementalscheme(scheme::AbstractScheme; kwargs...) -> Union{ElementalSchemes, AbstractChemical}

The elemental scheme of `scheme`.

A `Specific Methods` entry marked as `Entity` indicates that `elementalscheme` is applied for a scheme type that is not originally elemental.

# Generic Methods
* `AbstractScheme`: attribute seach.
* `AbstractStructuralScheme`: throw error.

# Specific Methods
* `ChemicalSchemes`: elemental scheme of all schemes.
* `IsotopomerizedSchemes`: elemental scheme of all schemes.

# Attribute Seach Workflow
1. Function: [`getchemicalattr`](@ref).
    * attribute: `:elementalscheme`.
    * default: itself.
"""
elementalscheme(scheme::AbstractScheme; kwargs...) = getchemicalattr(scheme, :elementalscheme, scheme, Union{ElementalSchemes, AbstractChemical})
elementalscheme(scheme::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("Structural scheme requires specific mapping to elemental scheme."))

"""
    structuralscheme(scheme::AbstractScheme; kwargs...) -> AbstractScheme

The structural scheme of `scheme`. All elemental schemes are regarded as structural schemes as well.

# Generic Methods
* `AbstractScheme`: attribute seach.

# Specific Methods
* `ChemicalSchemes`: structural scheme of all schemes.
* `IsotopomerizedSchemes`: structural scheme of all schemes.

# Attribute Seach Workflow
1. Function: [`getchemicalattr`](@ref).
    * attribute: `:structuralscheme`.
    * default: itself.
"""
structuralscheme(scheme::AbstractScheme; kwargs...) = getchemicalattr(scheme, :structuralscheme, scheme, AbstractScheme)

"""
    chemicalspecies(chemical::AbstractChemical; kwargs...) -> Vector{<:AbstractChemical}

The chemical species represented by a vector derived from a chemical entity or chemical species.

A `Specific Methods` entry marked as `Species` indicates that the return value may represent an attribute of a species.

# Generic Methods
* `AbstractChemical`: `[chemical]`.

# Specific Methods
* `Isobars`: `chemical.chemicals`.
* `ChemicalTransition{<:Isobars}`: a vector of chemical transition.
"""
chemicalspecies(chemical::AbstractChemical; kwargs...) = [chemical]

"""
    chemicaltransition(chemical::AbstractChemical; kwargs...) -> Vector{<:AbstractChemical}

The chemical entities that are analyzed in each stage of instrumental analysis.

A `Specific Methods` entry marked as `Transition` indicates that the return value may represent an attribute of a transition.

# Generic Methods
* `AbstractChemical`: `[chemical]`.

# Specific Methods
* `ChemicalTransition`: a vector of chemical entities.
* `Isobars`: the most abundant chemical transition.
"""
chemicaltransition(chemical::AbstractChemical; kwargs...) = [chemical]

"""
    chemicalparent(chemical::AbstractChemicalScheme; kwargs...) -> AbstractChemical

The parent chemical without delocalized isotopes replacement.

# Generic Methods
* `AbstractChemical`: itself.
* `AbstractScheme`: itself.
* `AbstractStructuralScheme`: throw error.

# Specific Methods
* Entity Level.
* `ChemicalTransition`: [`ChemicalTransition`](@ref) of parent chemicals of each transition.
* `Groupedisotopomers`: field `parent`.
* `ElementalScheme`: scheme of chemical parent of chemical entity.
* `ChemicalSchemes`: schemes of chemical parent(s) of chemical entity(s).
* `IsotopomerizedSchemes`: field `parent`.
* `Groupedisotopomerizedschemes`: field `parent`.
"""
chemicalparent(cc::AbstractChemical; kwargs...) = cc
chemicalparent(sch::AbstractScheme; kwargs...) = sch
chemicalparent(sch::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`chemicalparent` cannot be defined for structural scheme."))

"""
    isotopomersisotopes(chemical::AbstractChemicalScheme; loss = false, kwargs...) -> Vector{Pair{String, Int}}

The delocalized isotopes replacement of isotopomers.

# Generic Methods
* `AbstractChemical`: attribute seach.
* `AbstractScheme`: `Pair{String, Int}[]`.
* `AbstractStructuralScheme`: throw error.
* `AbstractCompleteScheme`: delocalized isotopes replacement of elemental scheme.

# Specific Methods
* Entity Level.
* `Isotopomers`: field `isotopes`.
* `Groupedisotopomers`: isotopes replacement of the most abundant isotopomers.
* `ElementalScheme`: change of delocalized isotopes replacement.
* `ChemicalSchemes`: all changes of delocalized isotopes replacement.
* `IsotopomerizedSchemes`: field `isotopes`.
* `Groupedisotopomerizedschemes`: isotopes replacement of the most abundant isotopomers.

# Attribute Seach Workflow
1. Function: [`getchemicalattr`](@ref).
    * attribute: `:isotopomersisotopes`.
    * default: `Pair{String, Int}[]`.

# Keyword arguments
* `loss::Bool` determines whether the chemical participates in a loss process, with sign flips propagated into the element counts.
"""
isotopomersisotopes(cc::AbstractChemicalScheme; loss = false, kwargs...) = reverse_elements(getchemicalattr(cc, :isotopomersisotopes, Pair{String, Int}[]), loss)
isotopomersisotopes(sch::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`isotopomersisotopes` cannot be defined for structural scheme."))
isotopomersisotopes(sch::AbstractCompleteScheme; loss = false, kwargs...) = isotopomersisotopes(elementalscheme(sch); loss, kwargs...) 

"""
    mass_shift_index(cc::AbstractChemical; isotope = "[13C]", ischemical = true, loss = false, kwargs...) -> Int 
    mass_shift_index(cc::AbstractScheme; isotope = "[13C]", ischemical = false, loss = false, kwargs...) -> Int 

A nominal index of mass shift between exact mass and monoisotopic mass using mass difference of `isotope` and its parent. 

# Generic Methods
* `AbstractChemicalScheme`: `mass_shift_index` calculated from `isotopomersisotopes`.

# Specific Methods
* Entity Level.

# Keyword arguments
* `ischemical::Bool` determines whether the input represents a chemical or a scheme.
* `loss::Bool` determines whether the chemical participates in a loss process, with sign flips propagated into `isotopomersisotopes`.
"""
mass_shift_index(cc::AbstractChemical; isotope_unit = nothing, isotope = "[13C]", loss = false, ischemical = true, kwargs...) = _mass_shift_index(isotopomersisotopes(cc), isnothing(isotope_unit) ? elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]] : isotope_unit; loss, ischemical, kwargs...)
mass_shift_index(cc::AbstractScheme; isotope_unit = nothing, isotope = "[13C]", loss = false, ischemical = false, kwargs...) = _mass_shift_index(isotopomersisotopes(cc), isnothing(isotope_unit) ? elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]] : isotope_unit; loss, ischemical, kwargs...)

@deprecate isotopomerstate mass_shift_index
"""
    groupedisotopomersisotopes(chemical::AbstractChemicalScheme; loss = false, kwargs...) -> Vector{Vector{Pair{String, Int}}}

The delocalized isotopes replacement of each isotopomers.

# Generic Methods
* `AbstractChemical`: single-element vector of `isotopomersisotopes`.
* `AbstractScheme`: `groupedisotopomersisotopes` of elemental scheme.
* `AbstractStructuralScheme`: throw error.

# Specific Methods
* Entity Level.
* `Groupedisotopomers`: isotopes replacements of each isotopomers.
* `ElementalScheme`: single-element vector of `isotopomersisotopes`.
* `ChemicalSchemes`: `Pair{String, Int}[]`.
* `IsotopomerizedSchemes`: single-element vector of `isotopomersisotopes`.
* `Groupedisotopomerizedschemes`: isotopes replacements of each isotopomers.

# Keyword arguments
* `loss::Bool` determines whether the chemical participates in a loss process, with sign flips propagated into the element counts.
"""
groupedisotopomersisotopes(x::AbstractChemical; loss = false, kwargs...) = [isotopomersisotopes(x; loss, kwargs...)]
groupedisotopomersisotopes(x::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`groupedisotopomersisotopes` cannot be defined for structural scheme."))
groupedisotopomersisotopes(x::AbstractScheme; loss = false, kwargs...) = groupedisotopomersisotopes(elementalscheme(x); loss, kwargs...)

"""
    groupedisotopomersabundance(chemical::AbstractChemicalScheme; loss = false, kwargs...) -> Vector{<:AbstractFloat}

The abundance of each isotopomers.

# Generic Methods
* `AbstractChemical`: `[1.0]`.
* `AbstractScheme`: `groupedisotopomersabundance` of elemental scheme.
* `AbstractStructuralScheme`: throw error.

# Specific Methods
* Entity Level.
* `Groupedisotopomers`: field `abundance`.
* `ElementalScheme`: `groupedisotopomersabundance` of field `chemical`.
* `ChemicalSchemes`: `[1.0]`.
* `IsotopomerizedSchemes`: `[1.0]`.
* `Groupedisotopomerizedschemes`: field `abundance`.
"""
groupedisotopomersabundance(x::AbstractChemical; kwargs...) = [1.0]
groupedisotopomersabundance(x::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`groupedisotopomersabundance` cannot be defined for structural scheme."))
groupedisotopomersabundance(x::AbstractScheme; kwargs...) = groupedisotopomersabundance(elementalscheme(x); kwargs...)
# truly formed chemical
"""
    analyzedchemical(chemical::AbstractChemicalScheme; kwargs...) -> AbstractChemical

The single chemical entity analyzed directly at the beginning of the instrumental analysis.

# Generic Methods
* `AbstractChemical`: itself.
* `AbstractScheme`: throw error.

# Specific Methods
* Species Level.
* `ChemicalTransition`: the very begining precursor, but cannot be chemical gain or loss.
* `Isobars`: the most abundant analyzed chemical.
"""
analyzedchemical(cc::AbstractChemical; kwargs...) = cc
analyzedchemical(cc::AbstractScheme; kwargs...) = throw(ArgumentError("Scheme cannot be analyzed directly."))

"""
    detectedchemical(chemical::AbstractChemicalScheme; kwargs...) -> AbstractChemical

The single chemical entity detected directly at the end of the instrumental analysis.

# Generic Methods
* `AbstractChemical`: itself or precursor with scheme (scheme is derived from `detectedchemical` of `precursor` and itself as product using [`completescheme`](@ref)).
* `AbstractScheme`: precursor with scheme. It requires a keyword argument `precursor`.

# Specific Methods
* Species Level.
* `ChemicalTransition`: the very ending product, but takes chemical gain and loss in consideration.
* `Isobars`: the most abundant detected chemical.
"""
function detectedchemical(cc::AbstractChemical; precursor = nothing, kwargs...)
    isnothing(precursor) && return cc 
    precursor = detectedchemical(precursor; kwargs...)
    detectedchemical(precursor, cc)
end
function detectedchemical(sch::AbstractScheme; precursor = nothing, kwargs...) 
    isnothing(precursor) && throw(ArgumentError("Scheme cannot be directly detected without precursor."))
    precursor = detectedchemical(precursor; kwargs...)
    detectedchemical(precursor, sch)
end
"""
    detectedisotopes(chemical::AbstractChemicalScheme; kwargs...) -> Vector{Pair{String, Int}}

The delocalized isotopes replacement of detected chemical. See [`detectedchemical`](@ref) for details.
"""
detectedisotopes(cc::AbstractChemicalScheme; kwargs...) = isotopomersisotopes(detectedchemical(cc); kwargs...)
function detectedisotopes(sch::AbstractScheme; precursor = nothing, precursorisotopes = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorisotopes) && throw(ArgumentError("Scheme cannot be directly detected without precursor information."))
    pre = isnothing(precursorisotopes) ? detectedisotopes(precursor; kwargs...) : precursorisotopes
    gain_elements(pre, isotopomersisotopes(sch; kwargs...))
end

"""
    detectedcharge(chemical::AbstractChemicalScheme; kwargs...) -> Int

The charge state of detected chemical. See [`detectedchemical`](@ref) for details.
"""
detectedcharge(cc::AbstractChemicalScheme; kwargs...) = charge(detectedchemical(cc); kwargs...)
function detectedcharge(sch::AbstractScheme; precursor = nothing, precursorcharge = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorcharge) && throw(ArgumentError("Scheme cannot be directly detected without precursor information."))
    (isnothing(precursorcharge) ? detectedcharge(precursor; kwargs...) : precursorcharge) + charge(sch; kwargs...)
end

# MS representation of chemical
"""
    inputchemical(chemical::AbstractChemicalScheme; kwargs...) -> AbstractChemical

The single chemical entity that serves as the input at the beginning of the instrumental analysis.

It is equivalent to [`analyzedchemical`](@ref), except that schemes are accepted. 
"""
inputchemical(cc::AbstractChemicalScheme; kwargs...) = cc

"""
    outputchemical(chemical::AbstractChemicalScheme; kwargs...) -> AbstractChemicalScheme

The single chemical entity that is the output at the end of the instrumental analysis.

It is equivalent to [`detectedchemical`](@ref), except that schemes are kept unchanged. 
"""
outputchemical(cc::AbstractChemicalScheme; kwargs...) = cc

"""
    seriesanalyzedchemical(chemical::AbstractChemicalScheme; kwargs...) -> Vector{<:AbstractChemical}

The chemical entities that are directly analyzed in each stage of instrumental analysis. 

It is equivalent to [`chemicaltransition`](@ref), except that schemes are transformed to detected chemicals.

# Generic Methods
* `AbstractChemicalScheme`: `[chemical]`.

# Specific Methods
* Transition Level.
* `ChemicalTransition`: a vector of analyzed chemicals. See [`analyzedchemical`](@ref) for details.
"""
seriesanalyzedchemical(cc::AbstractChemicalScheme; kwargs...) = [cc]

"""
    seriesanalyzedisotopes(chemical::AbstractChemicalScheme; kwargs...) -> Vector{Vector{Pair{String, Int}}}

The delocalized isotopes replacement of serially analyzed chemical. See [`seriesanalyzedchemical`](@ref) for details.
"""
seriesanalyzedisotopes(cc::AbstractChemicalScheme; kwargs...) = [isotopomersisotopes(c; kwargs...) for c in seriesanalyzedchemical(cc)]

"""
    seriesanalyzedcharge(chemical::AbstractChemicalScheme; kwargs...) -> Vector{Int}

The charge states of serially analyzed chemical. See [`seriesanalyzedchemical`](@ref) for details.
"""
seriesanalyzedcharge(cc::AbstractChemicalScheme; kwargs...) = [charge(c; kwargs...) for c in seriesanalyzedchemical(cc)]

"""
    msstage(chemical::AbstractChemical; kwargs...) -> Int

Number of stages of MS the chemical has been through.

# Generic Methods
* `AbstractChemical`: `1`.

# Specific Methods
* `ChemicalTransition`: `length(chemical.transition)`.
"""
msstage(cc::AbstractChemical; kwargs...) = 1

"""
    mmi(chemical::AbstractChemicalsScheme; loss = false, kwargs...) -> AbstractFloat

The monoisotopic mass of `chemical`.

# Generic Methods
* `AbstractChemical`: calculated from [`chemicalelements`](@ref) and [`charge`](@ref).
* `AbstractScheme`: calculated from [`chemicalelements`](@ref) and [`charge`](@ref); change of monoisotopic mass.
* `AbstractCompleteScheme`: change of monoisotopic mass of elemental scheme.

# Specific Methods
* Entity Level.
* `Isobars`: weighted mean of monoisotopic masses.
* `Groupedisotopomers`: weighted mean of monoisotopic masses.
* `ChemicalTransition`: monoisotopic mass of `analyzedchemical`.
* `Groupedisotopomerizedschemes`: weighted mean of change of monoisotopic masses.

# Keyword arguments
* `loss::Bool` determines whether the chemical participates in a loss process, with sign flips propagated through the charge.
"""
mmi(cc::AbstractChemicalScheme; loss = false) = mmi(chemicalelements(cc), charge(cc); loss)
mmi(sch::AbstractCompleteScheme; loss = false, kwargs...) = mmi(elementalscheme(sch); loss, kwargs...)

"""
    molarmass(chemical::AbstractChemicalScheme; loss = false, kwargs...) -> AbstractFloat

The molar mass of `chemical`.

# Generic Methods
* `AbstractChemical`: calculated from [`chemicalelements`](@ref) and [`charge`](@ref).
* `AbstractScheme`: calculated from [`chemicalelements`](@ref) and [`charge`](@ref); change of molar mass.
* `AbstractCompleteScheme`: change of molar mass of elemental scheme.

# Specific Methods
* Entity Level.
* `Isobars`: weighted mean of molar masses.
* `Groupedisotopomers`: weighted mean of molar masses.
* `ChemicalTransition`: molar mass of `analyzedchemical`.
* `Groupedisotopomerizedschemes`: weighted mean of change of molar masses.

# Keyword arguments
* `loss::Bool` determines whether the chemical participates in a loss process, with sign flips propagated through the charge.
"""
molarmass(cc::AbstractChemicalScheme; loss = false) = molarmass(chemicalelements(cc), charge(cc); loss)
molarmass(sch::AbstractCompleteScheme; loss = false, kwargs...) = molarmass(elementalscheme(sch); loss, kwargs...)

"""
    mz(chemical::AbstractChemical[, adduct]; loss = false, kwargs...) -> AbstractFloat
    mz(chemical::AbstractScheme; loss = false, kwargs...) -> AbstractFloat

The mass to charge ratio (m/z) of charged chemical or chemical with adduct. It is equivalent to `mmi(charged_chemical) / ncharge(charged_chemical)`.

# Generic Methods
* `AbstractChemical`: calculated from [`mmi`](@ref) and [`charge`](@ref).
* `AbstractAdductIon`: calculated from [`mmi`](@ref) and [`charge`](@ref) of both core chemical and adduct.
* `AbstractScheme`: calculated from [`mmi`](@ref) and [`charge`](@ref); change of m/z.
* `AbstractCompleteScheme`: change of m/z of elemental scheme.

# Specific Methods
* Entity Level.
* `Isobars`: weighted mean of m/z.
* `Groupedisotopomers`: weighted mean of m/z.
* `ChemicalTransition`: m/z of `analyzedchemical`.
* `Groupedisotopomerizedschemes`: weighted mean of change of m/z.

# Keyword arguments
* `loss:Bool` determines whether the chemical is part of chemical loss, and sign flips are factored out from elements.
"""
mz(charged_cc::AbstractChemicalScheme; loss = false, kwargs...) = charge(charged_cc) == 0 ? NaN : mmi(charged_cc; loss, kwargs...) / ncharge(charged_cc)
mz(cc::AbstractChemical, adduct; loss = false, kwargs...) = mz(ionize(cc; parse_adduct(adduct)...); loss, kwargs...)
function mz(adduct_ion::AbstractAdductIon; loss = false, kwargs...) 
    (ncore(adduct_ion) * mmi(ioncore(adduct_ion); loss, kwargs...) + mmi(ionadduct(adduct_ion); loss, kwargs...)) / ncharge(adduct_ion)
end
mz(sch::AbstractCompleteScheme; loss = false, kwargs...) = mz(elementalscheme(sch); loss, kwargs...)


# list all propertys
# add new propertys
# function chemicalpropertynames() end