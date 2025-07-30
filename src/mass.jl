"""
    mmi(formula::AbstractString, net_charge = 0)
    mmi(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}, net_charge = 0)
    mmi(chemical::AbstractChemical)
    mmi(chemicalpair::ChemicalPair)

Monoisotopic mass of a formula, collection of elements pairs or chemical.

This function calculates weighted average of monoisotopic masses for multiple chemicals. 
"""
function mmi(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}, net_charge = 0)
    # Vector of el => #el
    weight = 0.0u"g"
    for (el, n) in elements
        weight += ELEMENTS[:MASS][string(el)] * n
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
mmi(isobars::Isobars) = mean(mmi.(chemicalelements(isobars), charge.(isobars.chemicals)), weights(isobars.abundance))
mmi(isotopomers::Isotopomers) = mmi(chemicalelements(isotopomers), charge(isotopomers))
mmi(cp::ChemicalPair) = mmi(cp.precursor) => mmi(cp.product)

"""
    molarmass(formula::AbstractString, net_charge = 0)
    molarmass(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}, net_charge = 0)
    molarmass(chemical::AbstractChemical)
    molarmass(chemicalpair::ChemicalPair)

Molar mass of a formula, collection of elements pairs or chemical.

This function calculates weighted average of molar masses for multiple chemicals. 
"""
function molarmass(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}, net_charge = 0)
    # Vector of el => #el
    weight = 0.0u"g"
    for (el, n) in elements
        if haskey(ELEMENTS[:ISOTOPES], string(el))
            for i in ELEMENTS[:ISOTOPES][string(el)]
                weight += ELEMENTS[:MASS][i] * n * ELEMENTS[:ABUNDANCE][i]
            end
        else
            weight += ELEMENTS[:MASS][string(el)] * n
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
molarmass(isobars::Isobars) = mean(molarmass.(chemicalelements(isobars), charge.(isobars.chemicals)), weights(isobars.abundance))
molarmass(isotopomers::Isotopomers) = mmi(chemicalelements(isotopomers), charge(isotopomers))
molarmass(cp::ChemicalPair) = molarmass(cp.precursor) => molarmass(cp.product)

"""
    mz(charged_chemical::AbstractChemical)
    mz(chemical::AbstractChemical, adduct::AbstractAdduct)
    mz(chemical::AbstractChemical, adduct::AbstractString)
    mz(adduct_ion::AbstractAdductIon, adduct = ionadduct(adduct_ion))
    mz(isobars::Isobars) 
    mz(isobars::Isobars, adduct) 
    mz(chemicalpair::ChemicalPair)

Calculate m/z (monoisotopic mass over # charge) of `charged_chemical`, `adduct_ion`, `chemicalpair`, `isobars` or `chemical` with adduct `adduct`. It is equivalent to `mmi(cc) / ncharge(cc)`.
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
mz(cp::ChemicalPair) = mz(cp.precursor) => mz(cp.product)

# function mz(m::AbstractChemical, adduct::AbstractString)
#     ad = get(ADDUCT_NAME, string(adduct), nothing)
#     isnothing(ad) || return mz(m, ad)
#     adduct = replace(adduct, ADDUCT_ABBR...)
#     ion, charge = split(adduct, "]")
#     charge = split(charge, "+", keepempty = false)
#     charge = isempty(charge) ? 1 : begin
#         charge = split(only(charge), "-", keepempty = false)
#         isempty(charge) ? 1 : parse(Int, only(charge))
#     end
#     ion = replace(ion, "[" => "")
#     # #chemicals
#     nm, ion = split(ion, "M")
#     nm = isempty(nm) ? 1 : parse(Int, nm)

#     pos_adds = split(ion, "+", keepempty = false)
#     madd = mapreduce(+, pos_adds) do pos_add
#         # M-H / M+FA-H
#         mapreduce(mmi, -, split(pos_add, "-"))
#     end
#     (mmi(m) * nm + madd) / charge
# end