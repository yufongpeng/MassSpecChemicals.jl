"""
    parse_chemical([parser::AbstractChemicalParser,] name::AbstractString...; kwargs...) 
    parse_chemical([parser::AbstractChemicalParser,] chemical::AbstractChemical...; kwargs...) -> AbstractChemical
    parse_chemical([parser::AbstractChemicalParser,] chemicals::AbstractVector; kwargs...) 
    parse_chemical([parser::AbstractChemicalParser,] pair::Pair kwargs...) 

Parse chemical name and construct a chemical object using `parser`. The default parser is `ChemicalTransitionParser()`.
"""
parse_chemical(name...; kwargs...) = parse_chemical(ChemicalTransitionParser(), name...; kwargs...)
function parse_chemical(chemicalparser::ChemicalParser, name::AbstractString; kwargs...) 
    property = copy(chemicalparser.property)
    ks = first.(property)
    for (k, v) in kwargs
        i = findfirst(==(k), ks)
        if isnothing(i) && !isnothing(v)
            push!(property, k => v)
            push!(ks, k)
        elseif isnothing(i)
            continue
        else
            property[i] = k => v 
        end
    end
    filter!(x -> !isnothing(last(x)), property)
    i = findfirst(==(:formula), first.(property))
    if isnothing(i)
        i = findfirst(==(:elements), first.(property))
        isnothing(i) && throw(ArgumentError("No formula or elements are given"))
    end 
    fe = last(property[i])
    deleteat!(property, i)
    Chemical(name, chemicalelements(fe), property)
end

function parse_chemical(chemicalparser::FormulaChemicalParser, name::AbstractString; kwargs...) 
    property = copy(chemicalparser.property)
    ks = first.(property)
    for (k, v) in kwargs
        i = findfirst(==(k), ks)
        if isnothing(i) && !isnothing(v)
            push!(property, k => v)
            push!(ks, k)
        elseif isnothing(i)
            continue
        else
            property[i] = k => v 
        end
    end
    filter!(x -> !isnothing(last(x)), property)
    FormulaChemical(chemicalelements(name), property)
end

function parse_chemical(adductionparser::AdductIonParser, name::AbstractString; charge = nothing, kwargs...) 
    charge = isnothing(charge) ? adductionparser.charge : charge
    m = match(r"^\[(\d*)(.*)\](\d*[+-])*$", name)
    if isnothing(m) 
        n, f = match(r"^(\d*)(.*)$", name) 
        c = nothing
    else
        n, f, c = m 
        isnothing(c) && (charge = 0) 
    end
    n = isempty(n) ? 0 : parse(Int, n)
    if isnothing(c) && abs(charge) == 0 
        parse_chemical(adductionparser.chemicalparser, string(n, f); charge, kwargs...)
    elseif isnothing(c)
        nm = n > 1 ? n : ""
        lm = abs(charge) > 1 ? abs(charge) : ""
        a = charge > 0 ? string("[", nm, "M]", lm , "+") : string("[", nm, "M]", lm, "-")
        AdductIon(parse_chemical(adductionparser.chemicalparser, f; charge = nothing, kwargs...), parse_adduct(adductionparser.adductparser, a))
    else 
        plusa = split(f, "+"; limit = 2)
        minusa = split(popfirst!(plusa), "-"; limit = 2)
        m = popfirst!(minusa)
        a = string("[", n > 1 ? n : "", "M", isempty(minusa) ? "" : string("-", first(minusa)), isempty(plusa) ? "" : string("+", first(plusa)), "]", c)
        AdductIon(parse_chemical(adductionparser.chemicalparser, m; charge = nothing, kwargs...), parse_adduct(adductionparser.adductparser, a))
    end
end

function parse_chemical(chemicalparser::ChemicalGainLossParser, name::AbstractString; precursorcharge = nothing, kwargs...)
    precursorcharge = isnothing(precursorcharge) ? chemicalparser.charge : precursorcharge
    if startswith(name, "-") 
        loss = chemicalparser.loss
        if precursorcharge == 0 
            loss = 0
        else
            while (precursorcharge - loss) * precursorcharge <= 0 
                loss -= sign(precursorcharge)
            end
        end
        ChemicalLoss(parse_chemical(chemicalparser.chemicalparser, name[begin + 1:end]; charge = loss))
    elseif startswith(name, "+")
        gain = chemicalparser.gain
        if precursorcharge == 0 
            gain = 0
        else
            while (precursorcharge + gain) * precursorcharge <= 0 
                gain += sign(precursorcharge)
            end
        end
        ChemicalGain(parse_chemical(chemicalparser.chemicalparser, name[begin + 1:end]; charge = gain))
    else
        charge = min(abs(precursorcharge), abs(chemicalparser.charge)) * sign(precursorcharge)
        parse_chemical(chemicalparser.chemicalparser, name; charge)
    end
end

parse_chemical(::AbstractChemicalParser, cc::AbstractChemical...; kwargs...) = ChemicalSeries(cc...)

parse_chemical(::AbstractChemicalParser, cc::AbstractChemical, precursorcharge::Union{Int, Nothing}; kwargs...) = (cc, detectedcharge(cc; precursorcharge))
parse_chemical(::AbstractChemicalParser, cc::ChemicalTransition, precursorcharge::Union{Int, Nothing}; kwargs...) = (cc.transition, detectedcharge(cc))
function parse_chemical(chemicalparser::AbstractChemicalParser, name::AbstractString, precursorcharge::Union{Int, Nothing}; kwargs...) 
    chemical = parse_chemical(chemicalparser, name; precursorcharge, kwargs...)
    chemical, detectedcharge(chemical; precursorcharge)
end

parse_chemical(chemicalparser::ChemicalTransitionParser, name::AbstractString...; kwargs...) = parse_chemical(chemicalparser, collect(name); kwargs...)
parse_chemical(chemicalparser::ChemicalTransitionParser, name::AbstractString; kwargs...) = parse_chemical(chemicalparser, split(name, r"\s*->\s*"); kwargs...)
function parse_chemical(chemicalparser::ChemicalTransitionParser, name::AbstractVector; kwargs...) 
    trans = Any[]
    precursorcharge = nothing
    for nm in name
        c, precursorcharge = parse_chemical(chemicalparser.chemicalparser, nm, precursorcharge; kwargs...)
        push!(trans, c)
    end
    ChemicalSeries(trans...)
end

function parse_chemical(chemicalparser::ChemicalTransitionParser, name::Pair; kwargs...) 
    pre, precursorcharge = parse_chemical(chemicalparser, first(name), nothing; kwargs...)
    ChemicalSeries(first(pre_parse_chemical(chemicalparser, pre, last(name), precursorcharge; kwargs...)))
end

function pre_parse_chemical(chemicalparser::ChemicalTransitionParser, pre, name, precursorcharge; kwargs...)
    pre = vectorize(pre) 
    post, precursorcharge = parse_chemical(chemicalparser.chemicalparser, name, detectedcharge(last(pre); precursorcharge); kwargs...)
    vcat(pre, post), precursorcharge
end

function pre_parse_chemical(chemicalparser::ChemicalTransitionParser, pre, name::Pair, precursorcharge; kwargs...) 
    pre, precursorcharge = pre_parse_chemical(chemicalparser, pre, first(name), precursorcharge; kwargs...)
    pre_parse_chemical(chemicalparser, pre, last(name), precursorcharge; kwargs...)
end

function pre_parse_chemical(chemicalparser::ChemicalTransitionParser, pre, name::ChemicalTransition, precursorcharge; kwargs...) 
    post = Any[]
    for nm in chemicaltransition(name)
        c, precursorcharge = parse_chemical(chemicalparser.chemicalparser, nm, precursorcharge; kwargs...)
        push!(post, c)
    end
    vcat(pre, post...), precursorcharge
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
    print(io, "Co-eluting isobars: ", join([measure_name(first(x)) for x in ci.elution], " -> "), " -> ", join([msanalyzer_name(first(x)) for x in ci.msanalyzer], " -> "))
end