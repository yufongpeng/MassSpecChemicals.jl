decorator(::ElementalScheme{false}; loss = false, delim = "") = string(delim, loss ? "Gain_" : "Loss_")
decorator(::ElementalScheme{false, FormulaChemical}; loss = false, delim = "") = loss ? "+" : "-"
decorator(::ElementalScheme{false, <:Isotopomers{FormulaChemical}}; loss = false, delim = "") = loss ? "+" : "-"
decorator(::ElementalScheme{false, <:Isotopomers{<:AbstractAdductIon{<:FormulaChemical}}}; loss = false, delim = "") = loss ? "+" : "-"

decorator(::ElementalScheme{true}; loss = false, delim = "") = string(delim, loss ? "Loss_" : "Gain_")
decorator(::ElementalScheme{true, FormulaChemical}; loss = false, delim = "") = loss ? "-" : "+"
decorator(::ElementalScheme{true, <:Isotopomers{FormulaChemical}}; loss = false, delim = "") = loss ? "-" : "+"
decorator(::ElementalScheme{true, <:Isotopomers{<:AbstractAdductIon{<:FormulaChemical}}}; loss = false, delim = "") = loss ? "-" : "+"

charge_repr(c) = c == 0 ? "" : string(abs(c) > 1 ? abs(c) : "", c > 0 ? "+" : "-") 
isotope_repr(isotopes) = isempty(unique_elements(isotopes)) ? "" : string("[", replace(chemicalformula(isotopes; delim = ","), "[" => "", "]" => ""), "]")

chemicalname(cc::Chemical; n = 1, kwargs...) = n > 1 ? string(n, cc.name) : cc.name
function chemicalname(cc::FormulaChemical; n = 1, bracket = false, kwargs...) 
    nm = bracket ? string("[", chemicalformula(cc; kwargs...), "]", charge_repr(charge(cc))) : chemicalformula(cc; kwargs...)
    n > 1 ? string(n, nm) : nm 
end
chemicalname(isobars::Isobars; verbose = true, kwargs...) = (length(chemicalspecies(isobars)) == 1 || verbose) ? string("Isobars[", join(chemicalname.(chemicalspecies(isobars); kwargs...), ", "), "]") : string("Isobars[", chemicalname(first(chemicalspecies(isobars); kwargs...)), ", …]")
chemicalname(isotopomers::Isotopomers; n = 1, kwargs...) = string(chemicalname(chemicalparent(isotopomers); n, kwargs...), isotope_repr(isotopomers.isotopes))

chemicalname(isotopomers::Groupedisotopomers; n = 1, kwargs...) = string(chemicalname(chemicalparent(isotopomers); n, kwargs...), isotopomers.state > 0 ? string("(+", isotopomers.state, ")") : isotopomers.state < 0 ? string("(-", abs(isotopomers.state), ")") : "") 
chemicalname(ct::ChemicalTransition; kwargs...) = join(chemicalname.(chemicaltransition(ct); kwargs...), " -> ")

chemicalname(sch::StructuralElementalScheme; n = 1, kwargs...) = chemicalname(elementalscheme(sch); n, kwargs...)
chemicalname(sch::ElementalScheme{true}; n = 1, loss = false, delim = "", kwargs...) = string(decorator(sch; loss, delim), chemicalname(chemicalentity(sch); n, loss, kwargs...))
chemicalname(sch::ElementalScheme{false}; n = 1, loss = false, delim = "", kwargs...) = string(decorator(sch; loss, delim), chemicalname(chemicalentity(sch); n, loss = !loss, kwargs...))
chemicalname(sch::IsotopomerizedSchema; n = 1, loss = false, bracket = false, delim = "|", kwargs...) = string(chemicalname(chemicalparent(sch); n, loss, bracket, delim, kwargs...), isotope_repr(sch.isotopes))
function chemicalname(sch::ChemicalSchema; n = 1, loss = false, bracket = false, delim = "|", kwargs...) 
    s = join([chemicalname(k; n = n * v, loss, bracket = false, delim, kwargs...) for (k, v) in pairs(sch.schema)], "")
    if startswith(s, delim)
        s = s[begin + 1:end] 
    end
    bracket ? string("[", s, "]") : s 
end

chemicalabbr(isobars::Isobars; verbose = true, kwargs...) = (length(chemicalspecies(isobars)) == 1 || verbose) ? string("Isobars[", join(chemicalabbr.(chemicalspecies(isobars); kwargs...), ", "), "]") : string("Isobars[", chemicalabbr(first(chemicalspecies(isobars); kwargs...)), ", …]")
chemicalabbr(isotopomers::Isotopomers; n = 1, kwargs...) = string(chemicalabbr(chemicalparent(isotopomers); n, kwargs...), isempty(unique_elements(isotopomers.isotopes)) ? "" : string("[", chemicalformula(isotopomers.isotopes; delim = ","), "]"))
chemicalabbr(isotopomers::Groupedisotopomers; n = 1, kwargs...) = string(chemicalabbr(chemicalparent(isotopomers); n, kwargs...), isotopomers.state > 0 ? string("(+", isotopomers.state, ")") : isotopomers.state < 0 ? string("(-", abs(isotopomers.state), ")") : "") 
chemicalabbr(ct::ChemicalTransition; kwargs...) = join(chemicalabbr.(chemicaltransition(ct); kwargs...), " -> ")

chemicalabbr(sch::StructuralElementalScheme; n = 1, bracket = true, kwargs...) = chemicalabbr(elementalscheme(sch); n, bracket, kwargs...)
chemicalabbr(sch::ElementalScheme{true, <:Electron}; n = 1, bracket = true, loss = false, kwargs...) = ""
chemicalabbr(sch::ElementalScheme{false, <:Electron}; n = 1, bracket = true, loss = false, kwargs...) = ""
chemicalabbr(sch::ElementalScheme{true}; n = 1, bracket = true, loss = false, kwargs...) = string(loss ? "-" : "+", chemicalabbr(chemicalentity(sch); n, bracket, kwargs...))
chemicalabbr(sch::ElementalScheme{false}; n = 1, bracket = true, loss = false, kwargs...) = string(loss ? "+" : "-", chemicalabbr(chemicalentity(sch); n, bracket, loss = !loss, kwargs...))
chemicalabbr(sch::IsotopomerizedSchema; loss = false, bracket = true, n = 1, kwargs...) = string(chemicalabbr(chemicalparent(sch); n, loss, bracket, kwargs...), isotope_repr(sch.isotopes))
function chemicalabbr(sch::ChemicalSchema; loss = false, bracket = true, n = 1, kwargs...) 
    s = join([chemicalabbr(k; n = n * v, loss, bracket = false, kwargs...) for (k, v) in pairs(sch.schema)], "")
    bracket ? string("[", s, "]", charge_repr(-charge(sch))) : s 
end

chemicalsmiles(isobars::Isobars; kwargs...) = chemicalsmiles(chemicalentity(isobars); kwargs...)
chemicalsmiles(isotopomers::Isotopomers; kwargs...) = chemicalsmiles(chemicalparent(isotopomers); kwargs...)
chemicalsmiles(isotopomers::Groupedisotopomers; kwargs...) = chemicalsmiles(chemicalentity(isotopomers); kwargs...)
chemicalsmiles(ct::ChemicalTransition; kwargs...) = chemicalsmiles(chemicalentity(ct); kwargs...)

function Base.show(io::IO, cc::AbstractChemical)
    nm = chemicalname(cc)
    print(io, isempty(nm) ? chemicalformula(cc) : nm)
end

function Base.show(io::IO, sch::AbstractScheme) 
    print(io, chemicalname(sch))
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