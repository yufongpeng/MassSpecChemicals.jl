charge(isobars::Isobars; kwargs...) = charge(chemicalentity(isobars); kwargs...)
charge(isotopomers::Isotopomers; kwargs...) = charge(chemicalparent(isotopomers); kwargs...)
charge(isotopomers::Groupedisotopomers; kwargs...) = charge(chemicalparent(isotopomers); kwargs...)
charge(loss::ChemicalLoss; kwargs...) = charge(chemicalentity(loss); kwargs...)
charge(gain::ChemicalGain; kwargs...) = charge(chemicalentity(gain); kwargs...)
charge(ct::ChemicalTransition; kwargs...) = charge(chemicalentity(ct); kwargs...) 

retentiontime(isobars::Isobars; kwargs...) = _isobar_species_attr(retentiontime, isobars; kwargs...)
retentiontime(isotopomers::Isotopomers; kwargs...) = retentiontime(chemicalparent(isotopomers); kwargs...)
retentiontime(isotopomers::Groupedisotopomers; kwargs...) = retentiontime(chemicalparent(isotopomers); kwargs...)
retentiontime(loss::ChemicalLoss; kwargs...) = retentiontime(chemicalentity(loss); kwargs...)
retentiontime(gain::ChemicalGain; kwargs...) = retentiontime(chemicalentity(gain); kwargs...)
retentiontime(ct::ChemicalTransition; kwargs...) = retentiontime(chemicalentity(ct); kwargs...)

msstage(isobars::Isobars{<: ChemicalTransition}; kwargs...) = only(unique(msstage.(chemicalspecies(isobars); kwargs...)))
msstage(ct::ChemicalTransition; kwargs...) = length(ct.transition)

_isobar_species_attr(fn, isobars::Isobars, args...; kwargs...) = mean(fn.(chemicalspecies(isobars), args...), weights(isobars.abundance))
_isobar_species_attr(fn, isobars::Isobars{<: ChemicalTransition}, args...; kwargs...) = mean(fn.(chemicalspecies(isobars), args...), weights(isobars.abundance[begin:size(isobars.abundance, 1)]))

"""
    mmi(formula::AbstractString, net_charge = 0) -> AbstractFloat
    mmi(elements, net_charge = 0) -> AbstractFloat

Monoisotopic mass of formula, and elements.
"""
function mmi(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}, net_charge = 0)
    # Vector of el => #el
    weight = 0.0u"g"
    for (el, n) in elements
        weight += elements_mass()[string(el)] * n
    end
    weight -= net_charge * ME
    ustrip(weight)
end

mmi(elements::Dictionary, net_charge = 0) = mmi(pairs(elements), net_charge)

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
mmi(isotopomers::Groupedisotopomers) = mmi(chemicalparent(isotopomers)) + mean([sum([(mmi(x) - mmi(elements_parents()[x])) * n for (x, n) in iso]) for iso in isotopomers.isotopes], weights(isotopomers.abundance))
mmi(isotopomers::Isotopomers) = mmi(chemicalelements(isotopomers), charge(isotopomers))
mmi(ct::ChemicalTransition) = mmi(analyzedchemical(ct))

"""
    molarmass(formula::AbstractString, net_charge = 0) -> AbstractFloat
    molarmass(elements, net_charge = 0) -> AbstractFloat

Molar mass of formula, and elements.
"""
function molarmass(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}, net_charge = 0)
    # Vector of el => #el
    weight = 0.0u"g"
    for (el, n) in elements
        if haskey(elements_isotopes(), string(el))
            for i in elements_isotopes()[string(el)]
                weight += elements_mass()[i] * n * elements_abundunce()[i]
            end
        else
            weight += elements_mass()[string(el)] * n
        end
    end
    weight -= net_charge * ME
    ustrip(weight)
end

molarmass(elements::Dictionary, net_charge = 0) = molarmass(pairs(elements), net_charge)

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
molarmass(isotopomers::Groupedisotopomers) = molarmass(chemicalparent(isotopomers)) + mean([sum([(molarmass(x) - molarmass(elements_parents()[x])) * n for (x, n) in iso]) for iso in isotopomers.isotopes], weights(isotopomers.abundance))
molarmass(isotopomers::Isotopomers) = molarmass(chemicalelements(isotopomers), charge(isotopomers))
molarmass(ct::ChemicalTransition) = molarmass(analyzedchemical(ct))

mz(isotopomers::Isotopomers{<: AbstractAdductIon}, adduct) = mz(Isotopomers(AdductIon(ioncore(isotopomers.parent), parse_adduct(adduct)), isotopomers.isotopes))
mz(isobars::Isobars{<: AbstractAdductIon}) = _isobar_species_attr(mz, isobars)
mz(isobars::Isobars, adduct) = _isobar_species_attr(mz, isobars, adduct)
mz(isotopomers::Groupedisotopomers{<: AbstractAdductIon}) = mz(chemicalparent(isotopomers)) + mean([sum([(mmi(x) - mmi(elements_parents()[x])) * n for (x, n) in iso]) for iso in isotopomers.isotopes], weights(isotopomers.abundance))
mz(isotopomers::Groupedisotopomers, adduct) = mz(chemicalparent(isotopomers), adduct) + mean([sum([(mmi(x) - mmi(elements_parents()[x])) * n for (x, n) in iso]) for iso in isotopomers.isotopes], weights(isotopomers.abundance))
mz(ct::ChemicalTransition) = mz(analyzedchemical(ct))