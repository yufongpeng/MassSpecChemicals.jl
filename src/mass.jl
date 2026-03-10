"""
    mmi(formula::AbstractString, net_charge = 0)
    mmi(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}, net_charge = 0) -> AbstractFloat
    mmi(chemical::AbstractChemical)

Monoisotopic mass of formula, elements or chemical entity.

This function calculates weighted average of monoisotopic masses for `Isobars`. 
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

mmi(cc::AbstractChemical) = mmi(chemicalelements(cc), charge(cc))
mmi(isobars::Isobars) = mean([mmi(chemicalelements(x), charge(x)) for x in isobars.chemicals], weights(isobars.abundance))
mmi(isotopomers::Isotopomers) = mmi(chemicalelements(isotopomers), charge(isotopomers))
mmi(cp::ChemicalPair) = mmi(cp.precursor)

"""
    molarmass(formula::AbstractString, net_charge = 0)
    molarmass(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}, net_charge = 0) -> AbstractFloat
    molarmass(chemical::AbstractChemical)

Molar mass of formula, elements or chemical entity.

This function calculates weighted average of molar masses for `Isobars`. 
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

molarmass(cc::AbstractChemical) = molarmass(chemicalelements(cc), charge(cc))
molarmass(isobars::Isobars) = mean([molarmass(chemicalelements(x), charge(x)) for x in isobars.chemicals], weights(isobars.abundance))
molarmass(isotopomers::Isotopomers) = molarmass(chemicalelements(isotopomers), charge(isotopomers))
molarmass(cp::ChemicalPair) = molarmass(cp.precursor)

"""
    pair_molarmass(chemicalpair::ChemicalPair) -> Pair{<: AbstractFloat, <: AbstractFloat}

Molar mass of of precursor and product.
"""
pair_molarmass(cp::ChemicalPair) = molarmass(cp.precursor) => molarmass(cp.product)

"""
    mz(charged_chemical::AbstractChemical)
    mz(chemical::AbstractChemical, adduct::AbstractAdduct)
    mz(chemical::AbstractChemical, adduct::AbstractString)
    mz(adduct_ion::AbstractAdductIon, adduct = ionadduct(adduct_ion)) -> AbstractFloat
    mz(isobars::Isobars) 
    mz(isobars::Isobars, adduct)

Mass to charge ratio (m/z) of charged chemical or chemical with adduct. It is equivalent to `mmi(charged_chemical) / ncharge(charged_chemical)`.
"""
mz(charged_cc::AbstractChemical) = charge(charged_cc) == 0 ? NaN : mmi(charged_cc) / ncharge(charged_cc)
mz(cc::AbstractChemical, adduct) = mz(AdductIon(cc, parse_adduct(adduct)))
function mz(adduct_ion::AbstractAdductIon) 
    adduct = ionadduct(adduct_ion)
    (kmer(adduct) * mmi(ioncore(adduct_ion)) + mmi(adductelements(adduct_ion)) - charge(adduct_ion) * ustrip(ME)) / ncharge(adduct_ion)
end
mz(adduct_ion::AbstractAdductIon, adduct) = mz(AdductIon(ioncore(adduct_ion), parse_adduct(adduct)))
mz(isotopomers::Isotopomers{<: AbstractAdductIon}, adduct) = mz(Isotopomers(AdductIon(ioncore(isotopomers.parent), parse_adduct(adduct)), isotopomers.isotopes))
mz(isobars::Isobars{<: AbstractAdductIon}) = mean(mz.(isobars.chemicals), weights(isobars.abundance))
mz(isobars::Isobars, adduct) = mean(mz.(isobars.chemicals, adduct), weights(isobars.abundance))
mz(cp::ChemicalPair) = mz(cp.precursor)