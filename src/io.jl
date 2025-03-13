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

function Base.show(io::IO, ri::UnionInterval)
    print(io, ri.intervals[1])
    i = 1
    while i < length(ri.intervals)
        i += 1
        print(io, " ∪ ", ri.intervals[i])
    end
end

function Base.show(io::IO, ::MIME"text/plain", ri::UnionInterval)
    print(io, typeof(ri), ": ")
    print(io, ri)
end

function Base.show(io::IO, ri::RealInterval)
    @match ri.leftoperator begin
        (&<)    => print(io, "(", ri.lowerbound, ", ")
        (&<=)   => print(io, "[", ri.lowerbound, ", ")
        x       => print(io, ri.lowerbound, " ", x ," x ")
    end
    @match ri.rightoperator begin
        (&<)    => print(io, ri.upperbound, ")")
        (&<=)   => print(io, ri.upperbound, "]")
        x       => print(io, x, " ", ri.upperbound)
    end
end

function Base.show(io::IO, ::MIME"text/plain", ri::RealInterval)
    print(io, typeof(ri), ": ")
    print(io, ri)
end

Base.show(io::IO, ri::EmptyInterval) = 
    print(io, "∅")

Base.show(io::IO, ::MIME"text/plain", ri::EmptyInterval) = 
    print(io, typeof(ri), ": ()")