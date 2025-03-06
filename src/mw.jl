
"""
    mw(formula::AbstractString, net_charge = 0)
    mw(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}, net_charge = 0)
    mw(chemical::AbstractChemical)

Molecular weight of a formula, collection of elements pairs or chemical.

Molecular weight is referred to molecular or formula mass for single monoisotopic chemicals, and weighted average for multiple chemicals. 
"""
function mw(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}, net_charge = 0)
    # Vector of el => #el
    weight = 0.0u"g"
    for (el, n) in elements
        weight += ELEMENTS[:MW][string(el)] * n
    end
    weight -= net_charge * ME
    ustrip(weight)
end

mw(elements::Dictionary, net_charge = 0) = mw(pairs(elements), net_charge)

function mw(formula::AbstractString, net_charge = 0)
    # if any prefix number
    n = match(r"^[0-9]+", formula)
    if isnothing(n)
        mw(chemicalelements(formula), net_charge)
    else
        formula = replace(formula, n.match => "")
        mw(chemicalelements(formula), net_charge) * parse(Int, n.match)
    end
end

mw(cc::AbstractChemical) = mw(chemicalformula(cc), charge(cc))
mw(isobars::Isobars) = mean(mw.(chemicalformula(isobars), charge.(isobars.chemicals)), weights(isobars.abundance))
mw(isotopomers::Isotopomers) = mw(chemicalformula(isotopomers), charge(isotopomers))

"""
    mz(charged_chemical::AbstractChemical)
    mz(chemical::AbstractChemical, adduct::AbstractAdduct)
    mz(chemical::AbstractChemical, adduct::AbstractString)
    mz(adduct_ion::AbstractAdductIon, adduct = ionadduct(adduct_ion))
    mz(isobars::Isobars) 
    mz(isobars::Isobars, adduct) 

Calculate m/z of `charged_chemical`, `adduct_ion` `isobars` or `chemical` with adduct `adduct`. It is equivalent to `mw(cc) / ncharge(cc)`.
"""
mz(charged_cc::AbstractChemical) = charge(charged_cc) == 0 ? NaN : mw(charged_cc, charge(charged_cc)) / ncharge(charged_cc)
mz(cc::AbstractChemical, adduct) = mz(AdductIon(cc, parse_adduct(adduct)))
function mz(adduct_ion::AbstractAdductIon) 
    adduct = ionadduct(adduct_ion)
    (kmer(adduct) * mw(ioncore(adduct_ion)) + mw(adductelements(adduct_ion)) - charge(adduct_ion) * ustrip(ME)) / ncharge(adduct_ion)
end
mz(adduct_ion::AbstractAdductIon, adduct) = mz(AdductIon(ioncore(adduct_ion), parse_adduct(adduct)))
mz(isobars::Isobars{<: AbstractAdductIon}) = mean(mz.(isobars.chemicals), weights(isobars.abundance))
mz(isobars::Isobars, adduct) = mean(mz.(isobars.chemicals, adduct), weights(isobars.abundance))

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
#         mapreduce(mw, -, split(pos_add, "-"))
#     end
#     (mw(m) * nm + madd) / charge
# end