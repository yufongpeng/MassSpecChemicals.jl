
"""
    mw(formula::AbstractString)
    mw(chemical::AbstractChemical)

Molecular weight of a formula or chemical.
"""
function mw(elements::Vector)
    # Vector of el => #el
    weight = 0.0u"g"
    for (el, n) in elements
        weight += MW[string(el)] * n
    end
    ustrip(weight)
end
function mw(formula::AbstractString)
    # if any prefix number
    n = match(r"^[0-9]+", formula)
    if isnothing(n)
        mw(parse_compound(transform_isotope_repr(formula)))
    else
        formula = replace(formula, n.match => "")
        mw(parse_compound(transform_isotope_repr(formula))) * parse(Int, n.match)
    end
end

mw(m::AbstractChemical) = mw(chemicalformula(m))
mw(m::Isobars) = mean(mw.(chemicalformula(m)), weights(m.abundance))

"""
    mz(chemical::AbstractChemical, adduct::AbstractAdduct)
    mz(chemical::AbstractChemical, adduct::AbstractString)
    mz(ion::AbstractIon, adduct = ionadduct(ion))
    mz(chemical::Isobars) 

Calculate m/z of `ion` or `chemical` with adduct `adduct`.
"""
mz(m::AbstractChemical, adduct) = mz(Ion(m, parse_adduct(adduct)))
function mz(ion::AbstractIon) 
    adduct = ionadduct(ion)
    (kmer(adduct) * mw(ioncore(ion)) + mw(adductelement(ion))) / charge(adduct)
end
mz(ion::AbstractIon, adduct) = mz(Ion(ioncore(ion), parse_adduct(adduct)))
mz(m::Isobars{<: AbstractIon}) = mean(mz.(m.chemicals), weights(m.abundance))
mz(m::Isobars, adduct) = mean(mz.(m.chemicals, adduct), weights(m.abundance))

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