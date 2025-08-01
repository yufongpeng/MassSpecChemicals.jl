function Base.show(io::IO, cc::AbstractChemical)
    print(io, isempty(chemicalname(cc)) ? chemicalformula(cc) : chemicalname(cc))
end

function Base.show(io::IO, adduct::T) where {T <: AbstractAdduct}
    radd = adductformula(adduct)
    isnothing(radd) && return print(io, T)
    k = kmer(adduct) > 1 ? string("[", kmer(adduct), "M") : "[M"
    print(io, k, radd, "]", ncharge(adduct) > 1 ? ncharge(adduct) : "", adduct isa AbstractPosAdduct ? "+" : "-")
end

function Base.show(io::IO, adduct_ion::AbstractAdductIon)
    print(io, chemicalname(adduct_ion))
end

function Base.show(io::IO, isobars::Isobars)
    print(io, chemicalname(isobars))
end

function Base.show(io::IO, isotopomers::Isotopomers)
    print(io, chemicalname(isotopomers))
end

function Base.show(io::IO, loss::ChemicalLoss)
    print(io, chemicalname(loss))
end
function Base.show(io::IO, cp::ChemicalPair)
    print(io, chemicalname(cp))
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