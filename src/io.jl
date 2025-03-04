function Base.show(io::IO, m::AbstractChemical)
    print(io, isempty(chemicalname(m)) ? chemicalformula(m) : chemicalname(m))
end

function Base.show(io::IO, adduct::T) where {T <: AbstractAdduct}
    radd = adductformula(adduct)
    isnothing(radd) && return print(io, T)
    k = kmer(adduct) > 1 ? string("[", kmer(adduct), "M") : "[M"
    print(io, k, radd, "]", ncharge(adduct) > 1 ? ncharge(adduct) : "", adduct isa AbstractPosAdduct ? "+" : "-")
end

function Base.show(io::IO, ion::AbstractIon)
    print(io, chemicalname(ion))
end

function Base.show(io::IO, m::Isobars)
    print(io, chemicalname(m))
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
    print(io, typeof(ri), ":\n ")
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
    print(io, typeof(ri), ":\n ")
    print(io, ri)
end

Base.show(io::IO, ri::EmptyInterval) = 
    print(io, "∅")

Base.show(io::IO, ::MIME"text/plain", ri::EmptyInterval) = 
    print(io, typeof(ri), "()")