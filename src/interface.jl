# implement == and hash for concrete type 
==(x::Chemical, y::Chemical) = x.name == y.name && sort(x.elements) == sort(y.elements) && sort(x.property; by = first) == sort(y.property; by = first)
==(x::FormulaChemical, y::FormulaChemical) = sort(x.elements) == sort(y.elements) && sort(x.property; by = first) == sort(y.property; by = first)
==(x::ChemicalTransition, y::ChemicalTransition) = x.transition == y.transition
==(x::Isobars, y::Isobars) = x.chemicals == y.chemicals && all(splat(isapprox), zip(x.abundance, y.abundance))
==(x::Isotopomers, y::Isotopomers) = x.parent == y.parent && sort(x.isotopes) == sort(y.isotopes) 
==(x::Groupedisotopomers, y::Groupedisotopomers) = x.parent == y.parent && x.state == y.state && x.isotope == y.isotope && all(v -> ==(sort(first(v)), sort(last(v))), zip(sort(x.isotopes), sort(y.isotopes))) && all(splat(isapprox), zip(x.abundance, y.abundance))
==(x::AdductIon, y::AdductIon) = x.core == y.core && x.adduct == y.adduct && x.ncore == y.ncore
==(x::T, y::T) where {T <: AbstractChemicalWrapper} = x.chemical == y.chemical
==(x::ElementalScheme{T}, y::ElementalScheme{T}) where T = x.chemical == y.chemical
==(x::StructuralElementalScheme, y::StructuralElementalScheme) = x.structuralscheme == y.structuralscheme && x.elementalscheme == y.elementalscheme
==(x::ChemicalSchema, y::ChemicalSchema) = x.schema == y.schema 
==(x::IsotopomerizedSchema, y::IsotopomerizedSchema) = x.parent == y.parent && all(v -> ==(sort(first(v)), sort(last(v))), zip(sort(x.isotopes), sort(y.isotopes)))
==(x::Groupedisotopomerizedschema, y::Groupedisotopomerizedschema) = x.parent == y.parent && x.state == y.state && x.isotope == y.isotope && all(v -> ==(sort(first(v)), sort(last(v))), zip(sort(x.isotopes), sort(y.isotopes))) && all(splat(isapprox), zip(x.abundance, y.abundance))

const ABUNDANCE_ROUNDING_DIGITS = 4
abundance_rounding_digits() = ABUNDANCE_ROUNDING_DIGITS

function hash(x::Chemical, h::UInt) 
    h = hash(x.name, h) 
    h = hash(sum(hash, x.elements; init = zero(h)), h)
    h = hash(map(first, axes(x.property)), h)
    h = hash(map(last, axes(x.property)), h)
    hash(sum(hash, x.property; init = zero(h)), h)
end
function hash(x::FormulaChemical, h::UInt) 
    h = hash(sum(hash, x.elements; init = zero(h)), h)    
    h = hash(map(first, axes(x.property)), h)
    h = hash(map(last, axes(x.property)), h)
    hash(sum(hash, x.property; init = zero(h)), h)
end
hash(x::T, h::UInt) where {T<:ChemicalTransition} = hash(T, hash(x.transition, h))
function hash(x::Isobars, h::UInt) 
    h = hash(x.chemicals, h)
    h = hash(map(first, axes(x.abundance)), h)
    h = hash(map(last, axes(x.abundance)), h)
    for a in x.abundance
        h = hash(round(a; digits = abundance_rounding_digits()), h)
    end
    h
end
function hash(x::Isotopomers, h::UInt) 
    h = hash(x.parent, h) 
    h = hash(map(first, axes(x.isotopes)), h)
    h = hash(map(last, axes(x.isotopes)), h)
    hash(sum(hash, x.isotopes; init = zero(h)), h)
end
function hash(x::Groupedisotopomers, h::UInt) 
    h = hash(x.parent, h) 
    h = hash(x.state, h) 
    h = hash(x.isotope, h) 
    h = hash(map(first, axes(x.isotopes)), h)
    h = hash(map(last, axes(x.isotopes)), h)
    for y in x.isotopes 
        h = hash(map(first, axes(y)), h)
        h = hash(map(last, axes(y)), h)
        h = hash(sum(hash, y; init = zero(h)), h)
    end
    h = hash(map(first, axes(x.abundance)), h)
    h = hash(map(last, axes(x.abundance)), h)
    for a in x.abundance
        h = hash(round(a; digits = abundance_rounding_digits()), h)
    end
    h
end
function hash(x::AdductIon, h::UInt) 
    h = hash(x.core, h) 
    h = hash(x.adduct, h)
    hash(x.ncore, h)
end
hash(x::T, h::UInt) where {T<:AbstractChemicalWrapper} = hash(T, hash(x.chemical, h))
function hash(x::ElementalScheme{T}, h::UInt) where T
    h = hash(T, h) 
    hash(x.chemical, h)
end
function hash(x::StructuralElementalScheme, h::UInt) 
    h = hash(x.structuralscheme, h)
    hash(x.elementalscheme, h)
end
hash(x::T, h::UInt) where {T<:ChemicalSchema} = hash(T, hash(x.schema, h)) 
function hash(x::IsotopomerizedSchema, h::UInt) 
    h = hash(x.parent, h) 
    h = hash(map(first, axes(x.isotopes)), h)
    h = hash(map(last, axes(x.isotopes)), h)
    hash(sum(hash, x.isotopes; init = zero(h)), h)
end
function hash(x::Groupedisotopomerizedschema, h::UInt) 
    h = hash(x.parent, h) 
    h = hash(x.state, h) 
    h = hash(x.isotope, h) 
    h = hash(map(first, axes(x.isotopes)), h)
    h = hash(map(last, axes(x.isotopes)), h)
    for y in x.isotopes 
        h = hash(map(first, axes(y)), h)
        h = hash(map(last, axes(y)), h)
        h = hash(sum(hash, y; init = zero(h)), h)
    end
    h = hash(map(first, axes(x.abundance)), h)
    h = hash(map(last, axes(x.abundance)), h)
    for a in x.abundance
        h = hash(round(a; digits = abundance_rounding_digits()), h)
    end
    h
end
copy(x::Chemical) = Chemical(x.name, copy(x.elements), copy(x.property))
copy(x::FormulaChemical) = FormulaChemical(copy(x.elements), copy(x.property))
copy(x::ChemicalTransition) = ChemicalTransition(copy(x.transition))
copy(x::Isobars) = Isobars(copy(x.chemicals), copy(x.abundance))
copy(x::Isotopomers) = Isotopomers(copy(x.parent), copy(x.isotopes))
copy(x::Groupedisotopomers) = Groupedisotopomers(copy(x.parent), x.state, x.isotope, [copy(y) for y in x.isotopes], copy(x.abundance))
copy(x::AdductIon) = AdductIon(copy(x.core), copy(x.adduct), x.ncore) 
copy(x::T) where {T<:AbstractChemicalWrapper} = T(copy(x.chemical))
copy(x::ElementalScheme{T}) where T = ElementalScheme(T, copy(x.chemical))
copy(x::StructuralElementalScheme) = StructuralElementalScheme(copy(x.structuralscheme), copy(x.elementalscheme)) 
copy(x::ChemicalSchema) = ChemicalSchema(copy(x.schema)) 
copy(x::IsotopomerizedSchema) = IsotopomerizedSchema(copy(x.parent), copy(x.isotopes)) 
copy(x::Groupedisotopomerizedschema) = Groupedisotopomerizedschema(copy(x.parent), x.state, x.isotope, [copy(y) for y in x.isotopes], copy(x.abundance))

in(cc::AbstractChemical, isobars::Isobars) = any(i -> ischemicalequal(i, cc), isobars)
length(isobars::Isobars) = length(chemicalspecies(isobars))
length(cc::AbstractChemical) = 1
Broadcast.broadcastable(cc::AbstractChemical) = Ref(cc)

Broadcast.broadcastable(x::AbstractScheme) = Ref(x)

*(x::Criteria, y::Number) = Criteria(x.aval * y, x.rval * y)
/(x::Criteria, y::Number) = Criteria(x.aval / y, x.rval / y)

+(x::IntervalSet, y::Number) = IntervalSet([r + y for r in x.items])
-(x::IntervalSet, y::Number) = IntervalSet([r - y for r in x.items])
*(x::IntervalSet, y::Number) = IntervalSet([r * y for r in x.items])
/(x::IntervalSet, y::Number) = IntervalSet([r / y for r in x.items])
function +(x::T, y::Number) where {F, L <: Bound, R <: Bound, T <: Interval{F, L, R}}  
    if L == Unbounded
        f = nothing
    elseif isinf(x.first) && isinf(y) && y > 0
        f = x.first
    else
        f = x.first + y
    end
    if R == Unbounded
        l = nothing
    elseif isinf(x.last) && isinf(y) && y < 0
        l = x.last
    else
        l = x.last + y
    end
    LL = isnothing(f) ? Unbounded : isinf(f) ? Open : L
    RR = isnothing(l) ? Unbounded : isinf(l) ? Open : R
    Interval{F, LL, RR}(f, l)
end
function -(x::T, y::Number) where {F, L <: Bound, R <: Bound, T <: Interval{F, L, R}}  
    if L == Unbounded
        f = nothing
    elseif isinf(x.first) && isinf(y) && y < 0
        f = x.first
    else
        f = x.first - y
    end
    if R == Unbounded
        l = nothing
    elseif isinf(x.last) && isinf(y) && y > 0
        l = x.last
    else
        l = x.last - y
    end
    LL = isnothing(f) ? Unbounded : isinf(f) ? Open : L
    RR = isnothing(l) ? Unbounded : isinf(l) ? Open : R
    Interval{F, LL, RR}(f, l)
end
function *(x::T, y::Number) where {F, L, R, T <: Interval{F, L, R}}  
    if L == Unbounded
        f = nothing
    elseif isinf(x.first) && y == 0
        f = x.first
    else
        f = x.first * y
    end
    if R == Unbounded
        l = nothing
    elseif isinf(x.last) && y == 0
        l = x.last
    else
        l = x.last * y
    end
    LL = isnothing(f) ? Unbounded : isinf(f) ? Open : L
    RR = isnothing(l) ? Unbounded : isinf(l) ? Open : R
    if y < 0 
        Interval{F, RR, LL}(l, f)
    else
        Interval{F, LL, RR}(f, l)
    end
end
function /(x::T, y::Number) where {F, L , R, T <: Interval{F, L, R}} 
    if L == Unbounded
        f = nothing
    elseif isinf(x.first) && isinf(y)
        f = x.first
    else
        f = x.first / y
    end
    if R == Unbounded
        l = nothing
    elseif isinf(x.last) && isinf(y)
        l = x.last
    else
        l = x.last / y
    end
    LL = isnothing(f) ? Unbounded : isinf(f) ? Open : L
    RR = isnothing(l) ? Unbounded : isinf(l) ? Open : R
    if y < 0 
        Interval{F, RR, LL}(l, f)
    else
        Interval{F, LL, RR}(f, l)
    end
end