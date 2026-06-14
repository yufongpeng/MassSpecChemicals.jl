"""
    ischemicalequal(x::AbstractChemicalsSchema, y::AbstractChemicalsSchema) -> Bool

Determine whether two chemicals are chemically equivalent. By default, it transforms both chemicals by `ischemicalequaltransform` and compares them by `istransformedchemicalequal`.
"""
ischemicalequal(x::AbstractChemical, y::AbstractChemical) = istransformedchemicalequal(ischemicalequaltransform(x), ischemicalequaltransform(y))
ischemicalequal(x::AbstractScheme, y::AbstractScheme) = istransformedchemicalequal(ischemicalequaltransform(x), ischemicalequaltransform(y))
ischemicalequal(x::Isobars, y::Isobars) = istransformedchemicalequal(x, y)
ischemicalequal(x::Isotopomers, y::Isotopomers) = istransformedchemicalequal(x, y)
ischemicalequal(x::Groupedisotopomers, y::Groupedisotopomers) = istransformedchemicalequal(x, y)
ischemicalequal(x::ChemicalTransition, y::ChemicalTransition) = istransformedchemicalequal(x, y)
ischemicalequal(x::ChemicalSchema, y::ChemicalSchema) = istransformedchemicalequal(x, y)
ischemicalequal(x::StructuralElementalScheme, y::StructuralElementalScheme) = istransformedchemicalequal(x, y)
ischemicalequal(x::ElementalScheme{true}, y::ElementalScheme{true}) = istransformedchemicalequal(x, y)
ischemicalequal(x::ElementalScheme{false}, y::ElementalScheme{false}) = istransformedchemicalequal(x, y)
ischemicalequal(x::IsotopomerizedSchema, y::IsotopomerizedSchema) = istransformedchemicalequal(x, y)
ischemicalequal(x::Groupedisotopomerizedschema, y::Groupedisotopomerizedschema) = istransformedchemicalequal(x, y)
ischemicalequal(x::AbstractChemicalsSchema, y::AbstractChemicalsSchema) = istransformedchemicalequal(ischemicalequaltransform(x), ischemicalequaltransform(y))

"""
    ischemicalequaltransform(x::AbstractChemicalsSchema) -> AbstractChemicalsSchema

Return an object for comparison with other chemicals by `istransformedchemicalequal`. 
"""
ischemicalequaltransform(x::AbstractChemical) = x 
ischemicalequaltransform(x::AbstractScheme) = x 
ischemicalequaltransform(x::T) where {T <: AbstractChemicalWrapper} = ischemicalequaltransform(x.chemical)
ischemicalequaltransform(x::Isobars) = length(x) == 1 ? ischemicalequaltransform(chemicalentity(x)) : x
ischemicalequaltransform(x::Isotopomers) = isempty(unique_elements(x.isotopes)) ? x.parent : x 
ischemicalequaltransform(x::Groupedisotopomers) = length(x.isotopes) > 1 ? x : isempty(unique_elements(x.isotopes[1])) ? x.parent : Isotopomers(x.parent, x.isotopes[1]) 
ischemicalequaltransform(x::IsotopomerizedSchema) = isempty(unique_elements(x.isotopes)) ? x.parent : x 
ischemicalequaltransform(x::Groupedisotopomerizedschema) = length(x.isotopes) > 1 ? x : isempty(unique_elements(x.isotopes[1])) ? x.parent : Isotopomers(x.parent, x.isotopes[1]) 
ischemicalequaltransform(x::ElementalScheme{T}) where T = ElementalScheme(T, ischemicalequaltransform(x.chemical))
ischemicalequaltransform(x::StructuralElementalScheme) = StructuralElementalScheme(ischemicalequaltransform(x.structuralscheme), ischemicalequaltransform(x.elementalscheme))

"""
    istransformedchemicalequal(x::AbstractChemicalsSchema, y::AbstractChemicalsSchema) -> Bool

Determine whether two chemicals are chemically equivalent after applying `ischemicalequaltransform`. For `Chemical` and `FormulaChemical`, It tests the name and the elements composition. For others, it defaults to `isequal`.
"""
istransformedchemicalequal(x::AbstractChemicalsSchema, y::AbstractChemicalsSchema) = false
istransformedchemicalequal(x::AbstractChemical, y::AbstractChemical) = isequal(x, y)
istransformedchemicalequal(x::AbstractScheme, y::AbstractScheme) = isequal(x, y)
istransformedchemicalequal(x::AbstractAdductIon, y::AbstractAdductIon) = ischemicalequal(ionadduct(x), ionadduct(y)) && ischemicalequal(ioncore(x), ioncore(y))
istransformedchemicalequal(x::Chemical, y::Chemical) = 
    isequal(chemicalname(x), chemicalname(y)) && isequal(sort!(unique_elements(chemicalelements(x))), sort!(unique_elements(chemicalelements(y))))
istransformedchemicalequal(x::FormulaChemical, y::FormulaChemical) = 
    isequal(chemicalname(x), chemicalname(y)) 
istransformedchemicalequal(x::Isobars, y::Isobars) = all(ischemicalequal(a, b) for (a, b) in zip(x.chemicals, y.chemicals)) && all(isapprox(a, b) for (a, b) in zip(x.abundance, y.abundance))
istransformedchemicalequal(x::Isotopomers, y::Isotopomers) = ischemicalequal(x.parent, y.parent) && isequal(sort!(unique_elements(x.isotopes)), sort!(unique_elements(y.isotopes)))
istransformedchemicalequal(x::Groupedisotopomers, y::Groupedisotopomers) = ischemicalequal(x.parent, y.parent) && x.state == y.state && x.isotope == y.isotope && all(isequal(sort!(unique_elements(vx)), sort!(unique_elements(vy))) for (vx, vy) in zip(x.isotopes, y.isotopes)) && all(isapprox(a, b) for (a, b) in zip(x.abundance, y.abundance))
istransformedchemicalequal(x::ChemicalTransition, y::ChemicalTransition) = all(ischemicalequal.(x.transition, y.transition))

istransformedchemicalequal(x::Groupedisotopomerizedschema, y::Groupedisotopomerizedschema) = ischemicalequal(x.parent, y.parent) && x.state == y.state && x.isotope == y.isotope && all(isequal(sort!(unique_elements(vx)), sort!(unique_elements(vy))) for (vx, vy) in zip(x.isotopes, y.isotopes)) && all(isapprox(a, b) for (a, b) in zip(x.abundance, y.abundance))
istransformedchemicalequal(x::IsotopomerizedSchema, y::IsotopomerizedSchema) = istransformedchemicalequal(x.parent, y.parent) && isequal(sort!(unique_elements(x.isotopes)), sort!(unique_elements(y.isotopes)))
function istransformedchemicalequal(x::ChemicalSchema, y::ChemicalSchema) 
    uk = [false for _ in eachindex(y.schema)]
    for (kx, vx) in zip(x.schema, x.number)
        pass = false
        for (i, ky) in enumerate(y.schema)
            uk[i] && continue 
            if ischemicalequal(kx, ky) && vx == y.number[i]
                pass = true
                uk[i] = true
                continue
            end
        end
        pass || return false
    end
    true
end
istransformedchemicalequal(x::ElementalScheme, y::ElementalScheme) = ischemicalequal(x.chemical, y.chemical)
istransformedchemicalequal(x::StructuralElementalScheme, y::StructuralElementalScheme) = ischemicalequal(structuralscheme(x), structuralscheme(y)) && ischemicalequal(elementalscheme(x), elementalscheme(y))

"""
    ionize([constructor = AdductIon,] chemical; kwargs...) -> AbstractAdductIon

Ionize `chemical` and wrap with `constructor`.
"""
ionize(chemical::AbstractChemical; kwargs...) = ionize(AdductIon, chemical; kwargs...)
"""
    ionize(::Type{AdductIon}, chemical::AbstractChemical; adduct, ncore = 1, kwargs...) -> AdductIon

Ionize `chemical` and wrap with `AdductIon`.

# Keyword Arguments 
* `adduct::AbstractScheme`: adduct scheme.
* `ncore::Int`: number of core chemicals.
"""
ionize(::Type{AdductIon}, chemical::AbstractChemical; adduct, ncore = 1, kwargs...) = AdductIon(chemical, adduct, ncore)
ionize(::Type{AdductIon}, chemical::AbstractAdductIon; adduct, ncore = 1, kwargs...) = AdductIon(ioncore(chemical), adduct, ncore)
ionize(::Type{AdductIon}, chemical::Isobars; adduct, ncore = 1, kwargs...) = Isobars(ionize.(chemicalspecies(chemical); adduct, ncore, kwargs...), chemical.abundance)
ionize(::Type{AdductIon}, chemical::Isotopomers; adduct, ncore = 1, kwargs...) = Isotopomers(ionize(chemicalparent(chemical); adduct, ncore, kwargs...), chemical.isotopes)
ionize(::Type{AdductIon}, chemical::Groupedisotopomers; adduct, ncore = 1, kwargs...) = Groupedisotopomers(ionize(chemicalparent(chemical); adduct, ncore, kwargs...), chemical.state, chemical.isotope, chemical.isotopes, chemical.abundance)

"""
    isotopomerize(chemical::AbstractChemicalsSchema, isotopes) -> AbstractChemicalsSchema

Add delocalized isotopic replacements `isotopes` to `chemical`.
"""
isotopomerize(chemical::AbstractChemical, isotopes) = Isotopomers(chemical, isotopes)
isotopomerize(chemical::Isotopomers, isotopes) = Isotopomers(chemical.parent, collect(gain_elements(chemical.isotopes, isotopes)))
isotopomerize(sch::StructuralElementalScheme, isotopes) = StructuralElementalScheme(structuralscheme(sch), isotopomerize(elementalscheme(sch), isotopes))
isotopomerize(sch::ElementalScheme{true}, isotopes) = ElementalScheme(true, isotopomerize(sch.chemical, isotopes))
isotopomerize(sch::ElementalScheme{false}, isotopes) = ElementalScheme(false, isotopomerize(sch.chemical, reverse_elements(isotopes, true)))
isotopomerize(sch::ChemicalSchema, isotopes) = IsotopomerizedSchema(sch, isotopes)
isotopomerize(sch::IsotopomerizedSchema, isotopes) = IsotopomerizedSchema(sch.parent, collect(gain_elements(sch.isotopes, isotopes)))
isotopomerize(sch::T, isotopes) where {T<:AbstractScheme} = throw(ArgumentError("Cannot add isotopes information to $T."))

"""
    groupedisotopomerize(chemical::AbstractChemicalsSchema, isotopes) -> AbstractChemicalsSchema

Add delocalized isotopic replacements `isotopes` to `chemical`.
"""
groupedisotopomerize(chemical::AbstractChemical, state, isotope, isotopes, abundance) = Groupedisotopomers(chemical, state, isotope, isotopes, abundance)
groupedisotopomerize(chemical::Isotopomers, state, isotope, isotopes, abundance) = Groupedisotopomers(chemical.parent, state, isotope, isotopes, abundance)
groupedisotopomerize(chemical::Groupedisotopomers, state, isotope, isotopes, abundance) = Groupedisotopomers(chemical.parent, state, isotope, isotopes, abundance)
groupedisotopomerize(sch::StructuralElementalScheme, state, isotope, isotopes, abundance) = StructuralElementalScheme(structuralscheme(sch), groupedisotopomerize(elementalscheme(sch), state, isotope, isotopes, abundance))
groupedisotopomerize(sch::ElementalScheme{true}, state, isotope, isotopes, abundance) = ElementalScheme(true, groupedisotopomerize(sch.chemical, state, isotope, isotopes, abundance))
groupedisotopomerize(sch::ElementalScheme{false}, state, isotope, isotopes, abundance) = ElementalScheme(false, groupedisotopomerize(sch.chemical, state, isotope, reverse_elements.(isotopes, true), abundance))
groupedisotopomerize(sch::ChemicalSchema, state, isotope, isotopes, abundance) = Groupedisotopomerizedschema(sch, state, isotope, isotopes, abundance)
groupedisotopomerize(sch::IsotopomerizedSchema, state, isotope, isotopes, abundance) = Groupedisotopomerizedschema(sch.parent, state, isotope, isotopes, abundance)
groupedisotopomerize(sch::Groupedisotopomerizedschema, state, isotope, isotopes, abundance) = Groupedisotopomerizedschema(sch.parent, state, isotope, isotopes, abundance)
groupedisotopomerize(sch::T, state, isotope, isotopes, abundance) where {T<:AbstractScheme} = throw(ArgumentError("Cannot add isotopes information to $T."))

"""
    isgainscheme(sch::AbstractScheme) -> Bool

Whether `sch` contains only chemical gains.
"""
isgainscheme(sch) = false
isgainscheme(sch::ElementalScheme{true}) = true
isgainscheme(sch::ChemicalSchema) = all(isgainscheme, sch.schema)
isgainscheme(sch::IsotopomerizedSchema) = isgainscheme(sch.parent)

"""
    islossscheme(sch::AbstractScheme) -> Bool

Whether `sch` contains only chemical losses.
"""
islossscheme(sch) = false
islossscheme(sch::ElementalScheme{false}) = true
islossscheme(sch::ChemicalSchema) = all(islossscheme, sch.schema)
islossscheme(sch::IsotopomerizedSchema) = islossscheme(sch.parent) 

# completescheme(::chemicaltype, ::scheme)
# single precursor scheme dispatch: 
# completescheme(::AdductIon{chemicaltype, StructuralElementalScheme{structuraltype}}, ::scheme)
# Multiple precursor scheme dispatch: 
# completescheme(::AdductIon{chemicaltype, ChemicalSchema}, ::scheme)
"""
    completescheme(precursor::AbstractChemical, product::AbstractChemicalsSchema) -> CompleteSchema 

Transform `product` into `CompleteSchema` with `precursor`. 

For generic chemical types, this function calls `structure_search` which searches the property `:structure` of `ioncore(precursor)` for `ionadduct(precursor)` and `product`. 
For other chemicals, it returns a complete scheme directly without incorporating any information from `precursor`.

New method should be defined for new structural scheme and new chemical type to generate valid complete scheme.
"""
completescheme(precursor::AbstractChemical, product::AbstractChemical) = StructuralElementalScheme(RandomProductScheme(), product)
completescheme(precursor::T, product::S) where {T<:AbstractChemical, S<:AbstractScheme} = throw(ArgumentError("Specific `completescheme(precursor::$T, scheme::$S)` method has to be implemented."))
completescheme(precursor::AbstractChemical, product::CompleteSchema) = product
completescheme(precursor::AbstractChemical, product::ChemicalSchema) = ChemicalSchema(completescheme.(Ref(precursor), product.schema), product.number)
completescheme(precursor::AbstractChemical, product::IsotopomerizedSchema) = IsotopomerizedSchema(completescheme(precursor, product.parent), product.isotopes)
completescheme(precursor::AbstractChemical, product::Groupedisotopomerizedschema) = Groupedisotopomerizedschema(completescheme(precursor, product.parent), product.state, product.isotope, product.isotopes, product.abundance)
completescheme(precursor::AbstractChemical, product::AbstractElementalScheme) = StructuralElementalScheme(product, copy(product))
# completescheme(precursor::AbstractChemical, product::ElementalScheme{T}) where T = StructuralElementalScheme(product, ElementalScheme(T, product))
# structure search for generic types
completescheme(precursor::GenericChemical, product::GenericChemical) = structure_search(precursor, nothing, product)
completescheme(precursor::GenericChemical, product::AdductIon{<:GenericChemical}) = structure_search(precursor, nothing, product)
completescheme(precursor::GenericChemical, product::AbstractScheme) = structure_search(precursor, nothing, product)
completescheme(precursor::GenericChemical, product::CompleteSchema) = product
completescheme(precursor::GenericChemical, product::ChemicalSchema) = structure_search(precursor, nothing, product)
completescheme(precursor::GenericChemical, product::IsotopomerizedSchema) = structure_search(precursor, nothing, product)
completescheme(precursor::GenericChemical, product::Groupedisotopomerizedschema) = structure_search(precursor, nothing, product)
completescheme(precursor::GenericChemical, product::AbstractElementalScheme) = structure_search(precursor, nothing, product)
completescheme(precursor::AdductIon{<:GenericChemical}, product::GenericChemical) = structure_search(ioncore(precursor), ionadduct(precursor), product)
completescheme(precursor::AdductIon{<:GenericChemical}, product::AdductIon{<:GenericChemical}) = structure_search(ioncore(precursor), ionadduct(precursor), product)
completescheme(precursor::AdductIon{<:GenericChemical}, product::AbstractScheme) = structure_search(ioncore(precursor), ionadduct(precursor), product)
completescheme(precursor::AdductIon{<:GenericChemical}, product::CompleteSchema) = product
completescheme(precursor::AdductIon{<:GenericChemical}, product::IsotopomerizedSchema) = structure_search(ioncore(precursor), ionadduct(precursor), product)
completescheme(precursor::AdductIon{<:GenericChemical}, product::ChemicalSchema) = structure_search(ioncore(precursor), ionadduct(precursor), product)
completescheme(precursor::AdductIon{<:GenericChemical}, product::AbstractElementalScheme) = structure_search(ioncore(precursor), ionadduct(precursor), product)

# single precursor scheme dispatch: 
# adductionscheme(::AdductIon{chemicaltype, StructuralElementalScheme{structuraltype}}, ::StructuralElementalScheme{structuraltype})
# Multiple precursor scheme dispatch: 
# adductionscheme(::AdductIon{chemicaltype, ChemicalSchema}, ::StructuralElementalScheme{structuraltype})
# Multiple product scheme dispatch: 
# adductionscheme(::AdductIon{chemicaltype, StructuralElementalScheme{structuraltype}}, ::ChemicalSchema)
# Multiple precursor/product scheme dispatch: 
# adductionscheme(::AdductIon{chemicaltype, ChemicalSchema}, ::ChemicalSchema)
"""
    adductionscheme(precursor::AdductIon, product_scheme::CompleteSchema) -> CompleteSchema

Return a new scheme blending `ionadduct(precursor)` and `product_scheme`. 

For generic chemical types, this function calls `structure_search` which searches the property `:schema` of `ioncore(precursor)` for `ionadduct(precursor)` and `product`. 
For other chemicals, it returns a complete scheme directly without incorporating any information from `precursor`.

Defining new method is optional for complete scheme of new structural scheme and new chemical type unless scheme blending is intended to be implemented.
"""
adductionscheme(precursor::AdductIon, product::CompleteSchema) = ChemicalSchema(ionadduct(precursor), product)
# schema search for generic adduction
adductionscheme(precursor::AdductIon{<:GenericChemical}, product::CompleteSchema) = schema_search(ioncore(precursor), ionadduct(precursor), product)

# For customized adduction type of specific chemical type, implement 
# detectedchemical(::chemicaltype, ::CompleteSchema) -> adductiontype
# detectedchemical(::adductiontype, ::CompleteSchema) -> adductiontype
# adductionscheme(::adductiontype, ::CompleteSchema) -> CompleteSchema
# completescheme(::chemicaltype, ::AbstractScheme) -> CompleteSchema
# completescheme(::adductiontype, ::AbstractScheme) -> CompleteSchema

# Internal interfaces for detectedchemical
"""
    detectedchemical(precursor::AbstractChemical, product::AbstractChemicalsSchema) -> AbstractChemical 

The chemical directly detected in MS. Product of chemical entity is directly returned; `completescheme` is applied to scheme first. 
By default, `AdductIon` is returned for product generating from scheme, and any chemicals are directly returned. For `AdductIon`, `adductionscheme` is called for blending precursor and product scheme. 

Defining new method is optional unless for using other `AbstractAdductIon` type.
"""
detectedchemical(precursor::AbstractChemical, product::AbstractChemicalsSchema) = detectedchemical(precursor, completescheme(precursor, product))

detectedchemical(precursor::AbstractChemical, product::AbstractCompleteScheme{T, <:AbstractChemical}) where T = elementalscheme(product)
detectedchemical(precursor::AbstractChemical, product::CompleteSchema) = AdductIon(precursor, product, 1) 
detectedchemical(precursor::AbstractAdductIon, product::AbstractCompleteScheme{T, <:AbstractChemical}) where T = elementalscheme(product)
detectedchemical(precursor::T, product::CompleteSchema) where {T<:AbstractAdductIon} = throw(ArgumentError("Specific `detectedchemical(precursor::$T, scheme::StructuralElementalScheme)` method has to be implemented."))
detectedchemical(precursor::AdductIon, product::AbstractCompleteScheme{T, <:AbstractChemical}) where T = elementalscheme(product)
detectedchemical(precursor::AdductIon, product::CompleteSchema) = AdductIon(ioncore(precursor), adductionscheme(precursor, product), ncore(precursor))

detectedchemical(precursor::Isotopomers, product::AbstractChemical) = Isotopomers(product, Pair{String, Int}[])
detectedchemical(precursor::Isotopomers, product::Isotopomers) = product
function detectedchemical(precursor::Isotopomers, product::AbstractScheme) 
    chemical = detectedchemical(chemicalparent(precursor), chemicalparent(product))
    isotopes = gain_elements(isotopomersisotopes(precursor), isotopomersisotopes(product))
    Isotopomers(chemical, isotopes)
end
detectedchemical(precursor::Groupedisotopomers, product::AbstractChemical) = Groupedisotopomers(product, 0, precursor.isotope, groupedisotopomersisotopes(product), groupedisotopomersabundance(product))
detectedchemical(precursor::Groupedisotopomers, product::Isotopomers) = Groupedisotopomers(chemicalparent(product), isotopomerstate(product; isotope = precursor.isotope), precursor.isotope, groupedisotopomersisotopes(product), groupedisotopomersabundance(product))
detectedchemical(precursor::Groupedisotopomers, product::Groupedisotopomers) = product
function detectedchemical(precursor::Groupedisotopomers, product::AbstractScheme) 
    chemical = detectedchemical(chemicalparent(precursor), chemicalparent(product))
    isotopes = map((x, y) -> gain_elements(x, y), groupedisotopomersisotopes(precursor), groupedisotopomersisotopes(product))
    state = _isotopomerstate(first(isotopes), elements_mass()[precursor.isotope] - elements_mass()[elements_parents()[precursor.isotope]])
    Groupedisotopomers(chemical, state, precursor.isotope, isotopes, groupedisotopomersabundance(product))
end

chemicalentity(isobars::Isobars; kwargs...) = chemicalentity(first(chemicalspecies(isobars)))
chemicalentity(isotopomers::Groupedisotopomers; kwargs...) = Isotopomers(chemicalparent(isotopomers), isotopomersisotopes(isotopomers))
chemicalentity(ct::ChemicalTransition; kwargs...) = chemicalentity(first(chemicaltransition(ct)))

elementalscheme(sch::Groupedisotopomerizedschema; kwargs...) = Groupedisotopomerizedschema(elementalscheme(sch.parent; kwargs...), sch.state, sch.isotope, sch.isotopes, sch.abundance)
elementalscheme(sch::IsotopomerizedSchema; kwargs...) = IsotopomerizedSchema(elementalscheme(sch.parent; kwargs...), sch.isotopes)
elementalscheme(sch::ChemicalSchema; kwargs...) = ChemicalSchema(elementalscheme.(sch.schema; kwargs...), sch.number)
structuralscheme(sch::Groupedisotopomerizedschema; kwargs...) = Groupedisotopomerizedschema(structuralscheme(sch.parent; kwargs...), sch.state, sch.isotope, sch.isotopes, sch.abundance)
structuralscheme(sch::IsotopomerizedSchema; kwargs...) = IsotopomerizedSchema(structuralscheme(sch.parent; kwargs...), sch.isotopes)
structuralscheme(sch::ChemicalSchema; kwargs...) = ChemicalSchema(structuralscheme.(sch.schema; kwargs...), sch.number)
structuralscheme(::Nothing; kwargs...) = nothing 
elementalscheme(::Nothing; kwargs...) = nothing 

chemicalspecies(isobars::Isobars; kwargs...) = isobars.chemicals

function chemicaltransition(isobars::Isobars{<: ChemicalTransition}; kwargs...) 
    ct = chemicaltransition.(chemicalspecies(isobars))
    [Isobars(getindex.(ct, i), abundance) for (i, abundance) in enumerate(eachcol(isobars.abundance))]
end
chemicaltransition(ct::ChemicalTransition; kwargs...) = ct.transition

chemicalparent(isobars::Isobars; kwargs...) = chemicalparent(chemicalentity(isobars); kwargs...)
chemicalparent(isotopomers::Isotopomers; kwargs...) = isotopomers.parent 
chemicalparent(isotopomers::Groupedisotopomers; kwargs...) = isotopomers.parent 
chemicalparent(ct::ChemicalTransition; kwargs...) = ChemicalTransition(chemicalparent.(chemicaltransition(ct); kwargs...))

chemicalparent(sch::StructuralElementalScheme; kwargs...) = StructuralElementalScheme(structuralscheme(sch), chemicalparent(elementalscheme(sch); kwargs...))
chemicalparent(sch::ElementalScheme{T}; kwargs...) where T = ElementalScheme(T, chemicalparent(sch.chemical; kwargs...))
chemicalparent(sch::IsotopomerizedSchema; kwargs...) = sch.parent
chemicalparent(sch::ChemicalSchema; kwargs...) = ChemicalSchema(chemicalparent.(sch.schema; kwargs...), sch.number)
chemicalparent(sch::Groupedisotopomerizedschema; kwargs...) = sch.parent 

inputchemical(isobars::Isobars; kwargs...) = Isobars([inputchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, begin])
inputchemical(ct::ChemicalTransition; kwargs...) = first(chemicaltransition(ct))

outputchemical(isobars::Isobars; kwargs...) = Isobars([outputchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, begin])
outputchemical(ct::ChemicalTransition; kwargs...) = last(chemicaltransition(ct))

analyzedchemical(isobars::Isobars; kwargs...) = detectedchemical(isobars; kwargs...)
analyzedchemical(isobars::Isobars{<: ChemicalTransition}; kwargs...) = Isobars([analyzedchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, begin])
analyzedchemical(ct::ChemicalTransition; kwargs...) = detectedchemical(inputchemical(ct); kwargs...)

function seriesanalyzedchemical(isobars::Isobars{<: ChemicalTransition}; kwargs...) 
    ct = chemicaltransition.(chemicalspecies(isobars))
    [analyzedchemical(Isobars(getindex.(ct, i), abundance); kwargs...) for (i, abundance) in enumerate(eachcol(isobars.abundance))]
end
function seriesanalyzedchemical(ct::ChemicalTransition; kwargs...) 
    v = AbstractChemical[]
    precursor = nothing
    for c in chemicaltransition(ct) 
        push!(v, detectedchemical(c; precursor))
        precursor = last(v)
    end
    v
end
function seriesanalyzedisotopes(ct::ChemicalTransition; kwargs...)
    v = Vector{Pair{String, Int}}[]
    precursorisotopes = nothing
    for c in chemicaltransition(ct) 
        push!(v, detectedisotopes(c; precursorisotopes))
        precursorisotopes = last(v)
    end
    v
end
function seriesanalyzedcharge(ct::ChemicalTransition; kwargs...)
    v = Int[]
    precursorcharge = nothing
    for c in chemicaltransition(ct) 
        push!(v, detectedcharge(c; precursorcharge))
        precursorcharge = last(v)
    end
    v
end
function seriesanalyzedelements(ct::ChemicalTransition; kwargs...)
    v = Vector{Pair{String, Int}}[]
    precursorelements = nothing
    for c in chemicaltransition(ct) 
        push!(v, detectedelements(c; precursorelements))
        precursorelements = last(v)
    end
    v
end

detectedchemical(isobars::Isobars; kwargs...) = Isobars([detectedchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance)
detectedchemical(isobars::Isobars{<: ChemicalTransition}; kwargs...) = Isobars([detectedchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, end])

function detectedisotopes(sch::AbstractScheme; precursor = nothing, precursorisotopes = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorisotopes) && throw(ArgumentError("Scheme cannot be directly detected without precursor information."))
    pre = isnothing(precursorisotopes) ? detectedisotopes(precursor; kwargs...) : precursorisotopes
    gain_elements(pre, isotopomersisotopes(sch; kwargs...))
end
detectedisotopes(sch::StructuralElementalScheme{T, <:AbstractChemical}; precursor = nothing, precursorisotopes = nothing, kwargs...) where T = isotopomersisotopes(sch; kwargs...)
function detectedcharge(sch::AbstractScheme; precursor = nothing, precursorcharge = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorcharge) && throw(ArgumentError("Scheme cannot be directly detected without precursor information."))
    (isnothing(precursorcharge) ? detectedcharge(precursor; kwargs...) : precursorcharge) + charge(sch; kwargs...)
end
detectedcharge(sch::StructuralElementalScheme{T, <:AbstractChemical}; precursor = nothing, precursorcharge = nothing, kwargs...) where T = charge(sch; kwargs...)
function detectedelements(sch::AbstractScheme; precursor = nothing, precursorelements = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorelements) && throw(ArgumentError("Scheme cannot be directly detected without precursor information."))
    gain_elements(isnothing(precursorelements) ? detectedelements(precursor; kwargs...) : precursorelements, chemicalelements(sch; kwargs...))
end
detectedelements(sch::StructuralElementalScheme{T, <:AbstractChemical}; precursor = nothing, precursorelements = nothing, kwargs...) where T = chemicalelements(sch; kwargs...)

detectedchemical(ct::ChemicalTransition; kwargs...) = 
    detectedchemical(outputchemical(ct); precursor = @view(chemicaltransition(ct)[begin:end - 1]), kwargs...) 
detectedchemical(ct::AbstractVector; kwargs...) = 
    detectedchemical(last(ct); precursor = @view(ct[begin:end - 1]), kwargs...) 
detectedisotopes(ct::ChemicalTransition; kwargs...) = 
    detectedisotopes(outputchemical(ct); precursor = @view(chemicaltransition(ct)[begin:end - 1]), kwargs...) 
detectedisotopes(ct::AbstractVector; kwargs...) = 
    detectedisotopes(last(ct); precursor = @view(ct[begin:end - 1]), kwargs...) 
detectedcharge(ct::ChemicalTransition; kwargs...) = 
    detectedcharge(outputchemical(ct); precursor = @view(chemicaltransition(ct)[begin:end - 1]), kwargs...) 
detectedcharge(ct::AbstractVector; kwargs...) = 
    detectedcharge(last(ct); precursor = @view(ct[begin:end - 1]), kwargs...) 
detectedelements(ct::ChemicalTransition; kwargs...) = 
    detectedelements(outputchemical(ct); precursor = @view(chemicaltransition(ct)[begin:end - 1]), kwargs...) 
detectedelements(ct::AbstractVector; kwargs...) = 
    detectedelements(last(ct); precursor = @view(ct[begin:end - 1]), kwargs...) 

# GenericChemical property search
# structural -> complete, elemental -> complete
structure_search(chemical, precursor_schema, product_schema::CompleteSchema) = product_schema
structure_search(chemical, precursor_schema, product_schema::IsotopomerizedSchema) = IsotopomerizedSchema(structure_search(chemical, precursor_schema, product_schema.parent), product_schema.isotopes)
structure_search(chemical, precursor_schema, product_schema::ChemicalSchema) = ChemicalSchema([structure_search(chemical, precursor_schema, k) for k in product_schema.schema], product_schema.number)
function structure_search(chemical, precursor_schema, product::GenericChemical) 
    product_schema = chemicalscheme_search(product, nothing) 
    isnothing(product_schema) && return StructuralElementalScheme(RandomProductScheme(), product)
    sch = _structure_search(chemical, precursor_schema, product_schema)
    isnothing(sch) ? StructuralElementalScheme(RandomProductScheme(), product) : sch
end
function structure_search(chemical, precursor_schema, product::AdductIon{<:GenericChemical}) 
    product_schema = chemicalscheme_search(ioncore(product), ionadduct(product)) 
    isnothing(product_schema) && return StructuralElementalScheme(RandomProductScheme(), product)
    sch = _structure_search(chemical, precursor_schema, product_schema)
    isnothing(sch) ? StructuralElementalScheme(RandomProductScheme(), product) : sch
end
function structure_search(chemical, precursor_schema, product_schema::AbstractElementalScheme) 
    sch = _structure_search(chemical, precursor_schema, product_schema)
    isnothing(sch) ? StructuralElementalScheme(product_schema, copy(product_schema)) : sch
end

function _structure_search(chemical, precursor_schema, product_schema)
    schema = getchemicalproperty(chemical, :structure, nothing)
    isnothing(schema) && return nothing
    i = findfirst(x -> first(x) == structuralscheme(precursor_schema), schema)
    isnothing(i) && return nothing
    scheme = last(schema[i])
    i = findfirst(x -> first(x) == structuralscheme(product_schema), scheme)
    isnothing(i) && return nothing
    StructuralElementalScheme(first(scheme[i]), last(scheme[i]))
end

function chemicalscheme_search(chemical, scheme)
    schema = getchemicalproperty(chemical, :chemicalscheme, nothing)
    isnothing(schema) && return nothing
    i = findfirst(x -> first(x) == structuralscheme(scheme), schema)
    isnothing(i) && return nothing
    last(schema[i])
end

# complete -> complete
schema_search(chemical, precursor_schema::Nothing, product_schema::CompleteSchema) = product_schema
function schema_search(chemical, precursor_schema, product_schema::CompleteSchema)
    schema = getchemicalproperty(chemical, :schema, nothing)
    isnothing(schema) && return ChemicalSchema(precursor_schema, product_schema)
    i = findfirst(x -> first(x) == structuralscheme(precursor_schema), schema)
    isnothing(i) && return ChemicalSchema(precursor_schema, product_schema)
    scheme = last(schema[i])
    i = findfirst(x -> first(x) == structuralscheme(product_schema), scheme)
    isnothing(i) && return ChemicalSchema(precursor_schema, product_schema)
    structure_search(chemical, nothing, last(scheme[i]))
end

