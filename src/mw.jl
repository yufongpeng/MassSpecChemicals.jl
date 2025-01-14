
"""
    mw(formula::AbstractString)
    mw(m::AbstractChemical)

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
        mw(parse_compound(interpret_isotope(formula)))
    else
        formula = replace(formula, n.match => "")
        mw(parse_compound(interpret_isotope(formula))) * parse(Int, n.match)
    end
end

mw(m::AbstractChemical) = mw(chemicalformula(m))

"""
    mz(m::AbstractChemical, adduct::AbstractAdduct)
    mz(m::AbstractChemical, adduct::AbstractString)
    mz(ion::AbstractIon)

Calculate m/z of `ion` or `m` with adduct `adduct`.
"""
mz(m::AbstractChemical, adduct::AbstractAdduct) = adductmw(adduct)(mw(m))
mz(m::AbstractChemical, adduct::AbstractString) = adductmw(adduct)(mw(m))
mz(ion::AbstractIon) = adductmw(ion)(mw(ioncore(ion)))
mz(ions::IonCluster) = sum(ions.abundance .* mz.(ions.ions)) / sum(ions.abundance)

function adductmw(adduct::AbstractString)
    ad = get(DEFAULT_ADDUCT, string(adduct), nothing)
    isnothing(ad) || return adductmw(ad)
    adduct = replace(adduct, ADDUCT_FORMULA...)
    ion, charge = split(adduct, "]")
    charge = split(charge, "+", keepempty = false)
    charge = isempty(charge) ? 1 : begin
        charge = split(only(charge), "-", keepempty = false)
        isempty(charge) ? 1 : parse(Int, only(charge))
    end
    ion = replace(ion, "[" => "")
    # #chemicals
    nm, ion = split(ion, "M")
    nm = isempty(nm) ? 1 : parse(Int, nm)

    pos_adds = split(ion, "+", keepempty = false)
    madd = mapreduce(+, pos_adds) do pos_add
        # M-H / M+FA-H
        mapreduce(mw, -, split(pos_add, "-"))
    end
    x -> (x * nm + madd) / charge
end

adductmw(ion::AbstractIon) = adductmw(ionadduct(ion))
adductmw(adduct::AbstractAdduct) = (x -> (kmer(adduct) * x + mw(adductelement(adduct))) / charge(adduct))