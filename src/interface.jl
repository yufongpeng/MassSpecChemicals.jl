# implement copy, ==, and hash for concrete type 
_chemicalattributes_equal(x::Nothing, y::Nothing) = true
_chemicalattributes_equal(x::ChemicalAttribute, y::Nothing) = isnothing(x.attribute) || isempty(x.attribute)
_chemicalattributes_equal(x::Nothing, y::ChemicalAttribute) = isnothing(y.attribute) || isempty(y.attribute)
_chemicalattributes_equal(x::ChemicalAttribute, y::ChemicalAttribute) = x == y
function _chemicalelements_equal(x::Vector{Pair{String, Int}}, y::Vector{Pair{String, Int}}) 
    x = unique_elements(x)
    y = unique_elements(y)
    length(x) == length(y) || return false
    # all(in(y), x)
    sort!(x) == sort!(y)
end

function _chemicalelements_isless_equal(x::Vector{Pair{String, Int}}, y::Vector{Pair{String, Int}}) 
    x = unique_elements(x)
    y = unique_elements(y)
    if length(x) == length(y) 
        sort!(x) 
        sort!(y)
        x == y ? nothing : x < y
    else
        length(x) < length(y) 
    end
end

function ==(x::ChemicalAttribute, y::ChemicalAttribute)
    if isnothing(x.attribute) || isempty(x.attribute)
        return isnothing(y.attribute) || isempty(y.attribute)
    end 
    sort(x.attribute; by = first) == sort(y.attribute; by = first)
end

function ==(x::SchemeStructure, y::SchemeStructure)
    if isnothing(x.transformer) || isempty(x.transformer.vector)
        isnothing(y.transformer) || isempty(y.transformer.vector) || return false
    elseif isnothing(y.transformer) || isempty(y.transformer.vector)
        return false
    else
        x.transformer == y.transformer || return false
    end
    if isnothing(x.merger) || isempty(x.merger.vector)
        isnothing(y.merger) || isempty(y.merger.vector) || return false
    elseif isnothing(y.merger) || isempty(y.merger.vector)
        return false
    else
        x.merger == y.merger || return false
    end
    true
end

==(x::SchemeVector, y::SchemeVector) = x.vector == y.vector

isless(::Nothing, ::AbstractChemicalScheme) = true
isless(::AbstractChemicalScheme, ::Nothing) = false
function isless(x::AbstractChemicalScheme, y::AbstractChemicalScheme)
    nx = chemicalname(x) 
    ny = chemicalname(y)
    nx == ny || return nx < ny 
    leq = _chemicalelements_isless_equal(chemicalelements(x), chemicalelements(y)) 
    isnothing(leq) || return leq
    x == y ? false : hash_isless(x, y)
end

function hash_isless(x::S, y::T) where {S, T}
    px = propertynames(x)
    py = propertynames(y)
    ps = [px..., py...]
    ps = sort!(unique(ps))
    for p in ps
        pinx = p in px
        piny = p in py 
        if pinx && piny 
            hx = hash(getproperty(x, p))
            hy = hash(getproperty(y, p))
            hx == hy || return hx < hy
        else
            return piny
        end
    end
    false
end

==(x::Chemical, y::Chemical) =
    x.name == y.name && _chemicalelements_equal(x.elements, y.elements) && _chemicalattributes_equal(x.attribute, y.attribute)
==(x::FormulaChemical, y::FormulaChemical) =
    _chemicalelements_equal(x.elements, y.elements) && _chemicalattributes_equal(x.attribute, y.attribute)
function ==(x::T, y::T) where {T<:AttributedChemical}
    fields = fieldnames(T)
    attrname = getchemicalattrname(T)
    all(getfield(x, f) == getfield(y, f) for f in fields if f != :elements && f != attrname) &&
        _chemicalelements_equal(x.elements, y.elements) &&
        _chemicalattributes_equal(getchemicalattr(x, attrname, nothing), getchemicalattr(y, attrname, nothing))
end
==(x::ChemicalTransition, y::ChemicalTransition) = x.transition == y.transition
==(x::Isobars, y::Isobars) = x.chemicals == y.chemicals && all(splat(isapprox), zip(x.abundance, y.abundance))
==(x::Isotopomers, y::Isotopomers) = x.parent == y.parent && x.isotopes == y.isotopes
==(x::Groupedisotopomers, y::Groupedisotopomers) = x.parent == y.parent && x.index == y.index && x.isotope == y.isotope && all(splat(==), zip(x.isotopes, y.isotopes)) && all(splat(isapprox), zip(x.abundance, y.abundance))
==(x::AdductIon, y::AdductIon) = x.core == y.core && x.adduct == y.adduct && x.ncore == y.ncore
==(x::ElementalScheme{T}, y::ElementalScheme{T}) where T = x.chemical == y.chemical
==(x::StructuralElementalScheme, y::StructuralElementalScheme) = x.structuralscheme == y.structuralscheme && x.elementalscheme == y.elementalscheme
==(x::ChemicalSchemes, y::ChemicalSchemes) = x.schemes == y.schemes 
==(x::IsotopomerizedSchemes, y::IsotopomerizedSchemes) = x.parent == y.parent && x.isotopes == y.isotopes 
==(x::Groupedisotopomerizedschemes, y::Groupedisotopomerizedschemes) = x.parent == y.parent && x.index == y.index && x.isotope == y.isotope && all(splat(==), zip(x.isotopes, y.isotopes)) && all(splat(isapprox), zip(x.abundance, y.abundance))
function ==(x::ElementsVector, y::ElementsVector) 
    ixs = sortperm(x.elements)
    iys = sortperm(y.elements)
    ix = 1
    iy = 1
    @inbounds while true
        ix > lastindex(ixs) && iy > lastindex(iys) && break 
        if ix > lastindex(ixs)
            y.numbers[iys[iy]] == 0 || return false
            iy += 1
        elseif iy > lastindex(iys)
            x.numbers[ixs[ix]] == 0 || return false
            iy += 1
        elseif x.elements[ixs[ix]] == y.elements[iys[iy]]
            x.numbers[ixs[ix]] == y.numbers[iys[iy]] || return false
            ix += 1 
            iy += 1
        elseif x.numbers[ixs[ix]] == 0
            ix += 1
        elseif y.numbers[iys[iy]] == 0 
            iy += 1 
        else
            return false
        end
    end
    return true
end

const ABUNDANCE_ROUNDING_DIGITS = 4
abundance_rounding_digits() = ABUNDANCE_ROUNDING_DIGITS

function hash(x::SchemeVector, h::UInt)
    hash(x.vector, h) 
end
function hash(x::SchemeStructure, h::UInt)
    h = hash(SchemeStructure, h)
    h = hash(x.transformer, h) 
    hash(x.merger, h)
end
function hash(x::ChemicalAttribute, h::UInt)
    isnothing(x.attribute) ? hash(nothing, h) : hash(sort(x.attribute; by = first), h)
end
function hash(x::T, h::UInt) where {T<:AttributedChemical}
    h = hash(T, h)
    h = hash(sum(hash, sort_unique_elements(x.elements); init = zero(h)), h)
    attrname = getchemicalattrname(T)
    for f in fieldnames(T)
        if f != :elements && f != attrname
            h = hash(getfield(x, f), h)
        end
    end
    hash(getchemicalattr(x, attrname, nothing), h)
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
    hash(x.isotopes, h)
end
function hash(x::Groupedisotopomers, h::UInt) 
    h = hash(x.parent, h) 
    h = hash(x.index, h) 
    h = hash(x.isotope, h) 
    h = hash(map(first, axes(x.isotopes)), h)
    h = hash(map(last, axes(x.isotopes)), h)
    for y in x.isotopes
        h = hash(y, h)
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
function hash(x::ElementalScheme{T}, h::UInt) where T
    h = hash(T, h) 
    hash(x.chemical, h)
end
function hash(x::StructuralElementalScheme, h::UInt) 
    h = hash(x.structuralscheme, h)
    hash(x.elementalscheme, h)
end
hash(x::T, h::UInt) where {T<:ChemicalSchemes} = hash(T, hash(x.schemes, h)) 
function hash(x::IsotopomerizedSchemes, h::UInt) 
    h = hash(x.parent, h) 
    hash(x.isotopes, h)
end
function hash(x::Groupedisotopomerizedschemes, h::UInt) 
    h = hash(x.parent, h) 
    h = hash(x.index, h) 
    h = hash(x.isotope, h) 
    h = hash(map(first, axes(x.isotopes)), h)
    h = hash(map(last, axes(x.isotopes)), h)
    for y in x.isotopes 
        h = hash(y, h)
    end
    h = hash(map(first, axes(x.abundance)), h)
    h = hash(map(last, axes(x.abundance)), h)
    for a in x.abundance
        h = hash(round(a; digits = abundance_rounding_digits()), h)
    end
    h
end
function hash(x::ElementsVector, h::UInt)
    h = hash(map(first, axes(x.elements)), h)
    h = hash(map(last, axes(x.elements)), h)
    ixs = sortperm(x.elements)
    for ix in ixs
        h = hash(x.elements[ix], h)
        h = hash(x.numbers[ix], h)
    end
    h
end

stringcopy(x::T) where T = copy(x)
stringcopy(x::AbstractString) = x
stringcopy(x::Nothing) = x

copy(x::SchemeVector) = SchemeVector(stringcopy(x.vector))
copy(x::SchemeStructure) = SchemeStructure(stringcopy(x.transformer), stringcopy(x.merger))
copy(x::ChemicalAttribute) = ChemicalAttribute(stringcopy(x.attribute))
copy(x::T) where {T<:AbstractChemicalScheme} = T((stringcopy(getfield(x, f)) for f in fieldnames(T))...)
copy(x::ChemicalTransition) = ChemicalTransition(copy(x.transition))
copy(x::Isobars) = Isobars(copy(x.chemicals), copy(x.abundance))
copy(x::Isotopomers) = Isotopomers(copy(x.parent), copy(x.isotopes))
copy(x::Groupedisotopomers) = Groupedisotopomers(copy(x.parent), x.index, x.isotope, [copy(y) for y in x.isotopes], copy(x.abundance))
copy(x::AdductIon) = AdductIon(copy(x.core), copy(x.adduct), x.ncore) 
copy(x::ElementalScheme{T}) where T = ElementalScheme(T, copy(x.chemical))
copy(x::StructuralElementalScheme) = StructuralElementalScheme(copy(x.structuralscheme), copy(x.elementalscheme)) 
copy(x::ChemicalSchemes) = ChemicalSchemes(copy(x.schemes)) 
copy(x::IsotopomerizedSchemes) = IsotopomerizedSchemes(copy(x.parent), copy(x.isotopes)) 
copy(x::Groupedisotopomerizedschemes) = Groupedisotopomerizedschemes(copy(x.parent), x.index, x.isotope, [copy(y) for y in x.isotopes], copy(x.abundance))
copy(x::ElementsVector) = ElementsVector(copy(x.elements), copy(x.numbers))

iterate(ev::ElementsVector, i = 1) = length(ev.elements) < i ? nothing : (ev.elements[i] => ev.numbers[i], i + 1)
axes(ev::ElementsVector) = Base.OneTo(length(ev.elements))
isempty(ev::ElementsVector) = isempty(ev.elements)
length(ev::ElementsVector) = length(ev.elements)
eltype(ev::ElementsVector) = Pair{String, Int}

in(cc::AbstractChemical, isobars::Isobars) = any(i -> ischemicalequal(i, cc), isobars)
length(isobars::Isobars) = length(chemicalspecies(isobars))
length(cc::AbstractChemical) = 1
Broadcast.broadcastable(cc::AbstractChemical) = Ref(cc)
Broadcast.broadcastable(x::AbstractScheme) = Ref(x)

*(x::Criteria, y::Number) = Criteria(x.aval * y, x.rval * y)
/(x::Criteria, y::Number) = Criteria(x.aval / y, x.rval / y)
*(x::Var, y::Number) = Var(x.val * y)
/(x::Var, y::Number) = Var(x.val / y)
/(x::Var, y::Var) = x.val / y.val
iterate(x::Var, i = 1) = i == 1 ? x.val : nothing

+(x::IntervalSet, y::Number) = IntervalSet([r + y for r in x.items])
-(x::IntervalSet, y::Number) = IntervalSet([r - y for r in x.items])
*(x::IntervalSet, y::Number) = IntervalSet([r * y for r in x.items])
/(x::IntervalSet, y::Number) = IntervalSet([r / y for r in x.items])
function +(x::T, y::Number) where {F, L<:Bound, R<:Bound, T<:Interval{F, L, R}}  
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
function -(x::T, y::Number) where {F, L<:Bound, R<:Bound, T<:Interval{F, L, R}}  
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
function *(x::T, y::Number) where {F, L, R, T<:Interval{F, L, R}}  
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
function /(x::T, y::Number) where {F, L , R, T<:Interval{F, L, R}} 
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