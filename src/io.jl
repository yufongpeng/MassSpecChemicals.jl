"""
    parse_chemical(name; kwargs...)
    parse_chemical(::Type{Chemical}, name::AbstractString; kwargs...)

Parse chemical name and construct a chemical object. The default type is `Chemical`. `kwargs` are properties.
"""
parse_chemical(name; kwargs...) = parse_chemical(Chemical, name; kwargs...)
function parse_chemical(::Type{Chemical}, name::AbstractString; kwargs...) 
    ks = Dict{Symbol, Any}(kwargs)
    f = get!(ks, :elements, [])
    delete!(ks, :elements)
    if isempty(f) 
        f = get!(ks, :formula, "")
        delete!(ks, :formula)
        f = chemicalelements(f)
    end
    Chemical(name, f; ks...)
end

function parse_chemical(::Type{FormulaChemical}, name::AbstractString; charge = 0, kwargs...) 
    m = match(r"^\[(\d*)(.*)\](\d*[+-])*$", name)
    if isnothing(m) 
        n, f = match(r"^(\d*)(.*)$", name) 
        c = nothing
    else
        n, f, c = m 
    end
    n = isempty(n) ? 0 : parse(Int, n)
    if isnothing(c)
        if abs(charge) == 0 
            FormulaChemical(string(n, f); charge)
        else 
            nm = n > 1 ? n : ""
            lm = abs(charge) > 1 ? abs(charge) : ""
            a = charge > 0 ? string("[", nm, "M]", lm , "+") : string("[", nm, "M]", lm, "-")
            AdductIon(FormulaChemical(f), a)
        end
    else
        plusa = split(f, "+"; limit = 2)
        minusa = split(popfirst!(plusa), "-"; limit = 2)
        m = popfirst!(minusa)
        a = string("[", n > 1 ? n : "", "M", isempty(minusa) ? "" : string("-", first(minusa)), isempty(plusa) ? "" : string("+", first(plusa)), "]", c)
        AdductIon(FormulaChemical(m), a)
    end
end

decorator(loss::ChemicalLoss) = "Loss_"
decorator(loss::ChemicalLoss{FormulaChemical}) = "-"
decorator(loss::ChemicalLoss{<: Isotopomers{FormulaChemical}}) = "-"
decorator(loss::ChemicalLoss{<: Isotopomers{<: AbstractAdductIon{<: FormulaChemical}}}) = "-"

decorator(gain::ChemicalGain) = "Gain_"
decorator(gain::ChemicalGain{FormulaChemical}) = "+"
decorator(gain::ChemicalGain{<: Isotopomers{FormulaChemical}}) = "+"
decorator(gain::ChemicalGain{<: Isotopomers{<: AbstractAdductIon{<: FormulaChemical}}}) = "+"

chemicalname(cc::Chemical; kwargs...) = cc.name
chemicalname(cc::FormulaChemical; kwargs...) = chemicalformula(cc; kwargs...)
chemicalname(isobars::Isobars; verbose = true, kwargs...) = (length(chemicalspecies(isobars)) == 1 || verbose) ? string("Isobars[", join(chemicalname.(chemicalspecies(isobars); kwargs...), ", "), "]") : string("Isobars[", chemicalname(first(chemicalspecies(isobars); kwargs...)), ", …]")
chemicalname(isotopomers::Isotopomers; kwargs...) = string(chemicalname(chemicalparent(isotopomers); kwargs...), isempty(unique_elements(isotopomers.isotopes)) ? "" : string("[", replace(chemicalformula(isotopomers.isotopes; delim = ","), "[" => "", "]" => ""), "]"))
chemicalname(isotopomers::Isotopomers{FormulaChemical}; kwargs...) = string(chemicalformula(isotopomers))
function chemicalname(isotopomers::Isotopomers{<: AbstractAdductIon{<: FormulaChemical}}; kwargs...) 
    if isempty(adductelements(ionadduct(chemicalparent(isotopomers))))
        chemicalname(chemicalparent(isotopomers); corename = chemicalformula(isotopomers))
    else
        string(chemicalname(chemicalparent(isotopomers); kwargs...), isempty(unique_elements(isotopomers.isotopes)) ? "" : string("[", replace(chemicalformula(isotopomers.isotopes; delim = ","), "[" => "", "]" => ""), "]"))
    end
end
chemicalname(isotopomers::Groupedisotopomers; kwargs...) = string(chemicalname(chemicalparent(isotopomers)), isotopomers.state > 0 ? string("(+", isotopomers.state, ")") : isotopomers.state < 0 ? string("(-", abs(isotopomers.state), ")") : "") 
chemicalname(isotopomers::Groupedisotopomers{FormulaChemical}; kwargs...) = string(chemicalname(chemicalparent(isotopomers)), isotopomers.state > 0 ? string("(+", isotopomers.state, ")") : isotopomers.state < 0 ? string("(-", abs(isotopomers.state), ")") : "") 
chemicalname(loss::ChemicalLoss; kwargs...) = string(decorator(loss), chemicalname(chemicalentity(loss); kwargs...))
chemicalname(gain::ChemicalGain; kwargs...) = string(decorator(gain), chemicalname(chemicalentity(gain); kwargs...))
chemicalname(ct::ChemicalTransition; kwargs...) = join(chemicalname.(chemicaltransition(ct); kwargs...), " -> ")

chemicalabbr(isobars::Isobars; verbose = true, kwargs...) = (length(chemicalspecies(isobars)) == 1 || verbose) ? string("Isobars[", join(chemicalabbr.(chemicalspecies(isobars); kwargs...), ", "), "]") : string("Isobars[", chemicalabbr(first(chemicalspecies(isobars); kwargs...)), ", …]")
chemicalabbr(isotopomers::Isotopomers; kwargs...) = string(chemicalabbr(chemicalparent(isotopomers); kwargs...), isempty(unique_elements(isotopomers.isotopes)) ? "" : string("[", chemicalformula(isotopomers.isotopes; delim = ","), "]"))
chemicalabbr(isotopomers::Groupedisotopomers; kwargs...) = string(chemicalabbr(chemicalparent(isotopomers); kwargs...), isotopomers.state > 0 ? string("(+", isotopomers.state, ")") : isotopomers.state < 0 ? string("(-", abs(isotopomers.state), ")") : "") 
chemicalabbr(loss::ChemicalLoss; kwargs...) = string(decorator(loss), chemicalabbr(chemicalentity(loss); kwargs...))
chemicalabbr(gain::ChemicalGain; kwargs...) = string(decorator(gain), chemicalabbr(chemicalentity(gain); kwargs...))
chemicalabbr(ct::ChemicalTransition; kwargs...) = join(chemicalabbr.(chemicaltransition(ct); kwargs...), " -> ")

chemicalsmiles(isobars::Isobars; kwargs...) = chemicalsmiles(chemicalentity(isobars); kwargs...)
chemicalsmiles(isotopomers::Isotopomers; kwargs...) = chemicalsmiles(chemicalparent(isotopomers); kwargs...)
chemicalsmiles(isotopomers::Groupedisotopomers; kwargs...) = chemicalsmiles(chemicalentity(isotopomers); kwargs...)
chemicalsmiles(loss::ChemicalLoss; kwargs...) = chemicalsmiles(chemicalentity(loss); kwargs...) 
chemicalsmiles(gain::ChemicalGain; kwargs...) = chemicalsmiles(chemicalentity(gain); kwargs...) 
chemicalsmiles(ct::ChemicalTransition; kwargs...) = chemicalsmiles(chemicalentity(ct); kwargs...)

function Base.show(io::IO, cc::AbstractChemical)
    nm = chemicalname(cc)
    print(io, isempty(nm) ? chemicalformula(cc) : nm)
end

function Base.show(io::IO, adduct::T) where {T <: AbstractAdduct}
    radd = adductformula(adduct)
    isnothing(radd) && return print(io, T)
    k = kmer(adduct) > 1 ? string("[", kmer(adduct), "M") : "[M"
    print(io, k, radd, "]", ncharge(adduct) > 1 ? ncharge(adduct) : "", charge(adduct) > 0 ? "+" : "-")
end

repr_ri(ri::IntervalSet) = isempty(ri) ? "∅" : join([repr_ri(i) for i in ri.items], "∪")
function repr_ri(ri::Interval{T, L, R}) where {T, L, R}
    ff = L == Closed ? "[" : "("
    ll = R == Closed ? "]" : ")"
    f = L == Unbounded ? "-∞" : ri.first
    l = R == Unbounded ? "∞" : ri.last
    string(ff, f, ", ", l, ll)
end

function Base.show(io::IO, ri::IntervalSet)
    print(io, repr_ri(ri))
end

function Base.show(io::IO, c::Criteria{A, B}) where {A <: IntervalSet, B <: IntervalSet}
    print(io, "Criteria{IntervalSet, IntervalSet}(")
    print(io, repr_ri(c.aval), ", ")
    print(io, repr_ri(c.rval), ")")
end

function Base.show(io::IO, c::Criteria{A, B}) where {A <: Missing, B <: IntervalSet}
    print(io, "Criteria{Missing, IntervalSet}(")
    print(io, "missing, ")
    print(io, repr_ri(c.rval), ")")
end

function Base.show(io::IO, c::Criteria{A, B}) where {A <: IntervalSet, B <: Missing}
    print(io, "Criteria{IntervalSet, Missing}(")
    print(io, repr_ri(c.aval), ", ")
    print(io, "missing)")
end

function Base.show(io::IO, spec::Spectrum)
    print(io, "Spectrum of ", round(spec.initial_mass, digits = 4), " ~ ", round(spec.initial_mass + (length(spec.spectrum) - 1) * spec.binsize, digits = 4))
end

function Base.show(io::IO, ci::CoelutingIsobars)
    print(io, "Co-eluting isobars: ", join([colname(first(x)) for x in ci.elution], " -> "), " -> ", join([msanalyzername(first(x)) for x in ci.msanalyzer], " -> "))
end