charge(isobars::Isobars; kwargs...) = charge(chemicalentity(isobars); kwargs...)
charge(isotopomers::Isotopomers; kwargs...) = charge(chemicalparent(isotopomers); kwargs...)
charge(isotopomers::Groupedisotopomers; kwargs...) = charge(chemicalparent(isotopomers); kwargs...)
charge(ct::ChemicalTransition; kwargs...) = charge(chemicalentity(ct); kwargs...) 

charge(sch::AbstractCompleteScheme; kwargs...) = charge(elementalscheme(sch); kwargs...)
charge(loss::ElementalScheme{false}; kwargs...) = charge(loss.chemical; kwargs...)
charge(gain::ElementalScheme{true}; kwargs...) = -charge(gain.chemical; kwargs...)
charge(x::IsotopomerizedSchema) = charge(x.parent)
charge(x::ChemicalSchema) = sum(charge(k) * v for (k, v) in zip(x.schema, x.number))
charge(x::Groupedisotopomerizedschema) = charge(x.parent)

retentiontime(isobars::Isobars; kwargs...) = _isobar_species_attr(retentiontime, isobars; kwargs...)
retentiontime(isotopomers::Isotopomers; kwargs...) = retentiontime(chemicalparent(isotopomers); kwargs...)
retentiontime(isotopomers::Groupedisotopomers; kwargs...) = retentiontime(chemicalparent(isotopomers); kwargs...)
retentiontime(ct::ChemicalTransition; kwargs...) = retentiontime(chemicalentity(ct); kwargs...)

msstage(isobars::Isobars{<: ChemicalTransition}; kwargs...) = only(unique(msstage.(chemicalspecies(isobars); kwargs...)))
msstage(ct::ChemicalTransition; kwargs...) = length(ct.transition)

_isobar_species_attr(fn, isobars::Isobars, args...; kwargs...) = mean(fn.(chemicalspecies(isobars), args...), weights(isobars.abundance))

"""
    mmi(formula::AbstractString, net_charge = 0) -> AbstractFloat
    mmi(elements, net_charge = 0) -> AbstractFloat

Monoisotopic mass of formula, and elements.
"""
function mmi(elements::Union{<: Vector{<: Pair}, <: Dict}, net_charge = 0)
    # Vector of el => #el
    weight = 0.0
    for (el, n) in elements
        weight += elements_mass()[el] * n
    end
    weight -= net_charge * ME
    weight
end

function deltammi(elements::Union{<: Vector{<: Pair}, <: Dict}, net_charge = 0)
    weight = 0.0
    for (el, n) in elements
        weight += elements_mass()[el] * n
        weight -= elements_mass()[parent_element(el)] * n
    end
    weight -= net_charge * ME
    weight
end

function mmi(formula::AbstractString, net_charge = 0)
    # if any prefix number
    n = match(r"^[0-9]+", formula)
    if isnothing(n)
        mmi(chemicalelements(formula), net_charge)
    else
        formula = replace(formula, n.match => "")
        mmi(chemicalelements(formula), net_charge) * parse(Int, n.match)
    end
end

mmi(isobars::Isobars) = _isobar_species_attr(mmi, isobars)
mmi(isotopomers::Groupedisotopomers) = mmi(chemicalparent(isotopomers)) + mean([sum([(mmi(x) - mmi(elements_parents()[x])) * n for (x, n) in iso]; init = 0) for iso in isotopomers.isotopes], weights(isotopomers.abundance))
mmi(isotopomers::Isotopomers) = mmi(chemicalelements(isotopomers), charge(isotopomers))
mmi(ct::ChemicalTransition) = mmi(analyzedchemical(ct))

"""
    molarmass(formula::AbstractString, net_charge = 0) -> AbstractFloat
    molarmass(elements, net_charge = 0) -> AbstractFloat

Molar mass of formula, and elements.
"""
function molarmass(elements::Union{<: Vector{<: Pair}, <: Dict}, net_charge = 0)
    # Vector of el => #el
    weight = 0.0
    for (el, n) in elements
        if iselement(string(el))
            for i in elements_isotopes()[string(el)]
                weight += elements_mass()[i] * n * elements_abundance()[i]
            end
        else
            weight += elements_mass()[string(el)] * n
        end
    end
    weight -= net_charge * ME
    weight
end

function molarmass(formula::AbstractString, net_charge = 0)
    # if any prefix number
    n = match(r"^[0-9]+", formula)
    if isnothing(n)
        molarmass(chemicalelements(formula), net_charge)
    else
        formula = replace(formula, n.match => "")
        molarmass(chemicalelements(formula), net_charge) * parse(Int, n.match)
    end
end

molarmass(isobars::Isobars) = _isobar_species_attr(molarmass, isobars)
molarmass(isotopomers::Groupedisotopomers) = molarmass(chemicalparent(isotopomers)) + mean([sum([(molarmass(x) - molarmass(elements_parents()[x])) * n for (x, n) in iso]; init = 0) for iso in isotopomers.isotopes], weights(isotopomers.abundance))
molarmass(isotopomers::Isotopomers) = molarmass(chemicalelements(isotopomers), charge(isotopomers))
molarmass(ct::ChemicalTransition) = molarmass(analyzedchemical(ct))

mz(ct::ChemicalTransition) = mz(analyzedchemical(ct))
mz(ct::ChemicalTransition, adduct) = mz(analyzedchemical(ct), adduct)