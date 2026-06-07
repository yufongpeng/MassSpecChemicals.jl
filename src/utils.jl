"""
    match_chemical(exp, lib; colexp = :Chemical, collib = :Chemical) -> Table

Match chemicals in `exp` (a `Table` or `Vector`) to chemicals in `lib` (a `Table` or `Vector`). 
The resulting table is `exp` with matched index (column `LibID`), matched chemicals (column `Match`) and other information from `lib`.

If chemicals in `exp` are chemical pairs, the detected chemicals are utilized. See `detectedchemical` for details.
"""
function match_chemical(exp, lib; colexp = :Chemical, collib = :Chemical)
    del = Int[]
    libid = Int[]
    exp = hasproperty(exp, colexp) ? exp : Table(; Chemical = exp)
    chemical_exp = detectedchemical.(getproperty(exp, colexp))
    chemical_lib = hasproperty(lib, collib) ? getproperty(lib, collib) : lib
    for i in eachindex(exp)
        j = findfirst(x -> ischemicalequal(x, chemical_exp[i]), chemical_lib)
        isnothing(j) ? push!(del, i) : push!(libid, j)
    end
    id = setdiff(eachindex(exp), del)
    ps = filter(!=(collib), propertynames(lib))
    Table(exp[id]; LibID = libid, Match = chemical_lib[libid], [p => [getproperty(lib, p)[i] for i in libid] for p in ps]...)
end

"""
    acrit(x) -> Criteria

Cnstructing a `Criteria` with only absolute criteria.
"""
acrit(x) = Criteria(x, missing)
"""
    rcrit(x) -> Criteria

Cnstructing a `Criteria` with only relative criteria. 
"""
rcrit(x) = Criteria(missing, x)
"""
    crit(absolute)
    crit(x::Criteria)
    crit(absolute, relative) -> Criteria

Cnstructing a `Criteria`. When a non `Criteria` value is given, it consider it as absolute criteria; when two values are given, the first one is absolute, the second one is relative.
"""
crit(x) = acrit(x)
crit(x::Criteria) = x
crit(x, y) = Criteria(x, y)
# crit(x::Tuple) = Criteria(x...)

"""
    between(num::Number, range; lop = <=, rop = <=)
    between(num::Number; low, up, lop = <=, rop = <=)
    between(num::Number, value, tol; lop = <=, rop = <=)

Determine whether `num` lies in given range. 

The lower bound is `first(range)`, `low`, or `value - tol`; the upper bound is `last(range)`, `up`, or `value + tol`. 

`lop` and `rop` determine the operators for the lower bound and upper bound respectively.
"""
between(num::Number, range; lop = <=, rop = <=) = lop(first(range), num) && rop(num, last(range))
between(num::Number; low, up, lop = <=, rop = <=) = lop(low, num) && rop(num, up)
between(num::Number, value, tol; lop = <=, rop = <=) = lop(value - tol, num) && rop(num, value + tol)

"""
    @ri_str -> IntervalSet{Interval{Float64, L, R}}

Create a `IntervalSet{Interval{Float64, L, R}}` by mathematical real interval notation.

# Examples
```julia
julia> r = ri"[4, 7)"
1-interval IntervalSet{Interval{Float64, Closed, Open}}:
[4.0 .. 7.0)

julia> 4 in r
true

julia> 7 in r
false
```
"""
macro ri_str(expr)
    return @match expr begin 
        "()"    => intersect(IntervalSet(1.0..1.0), IntervalSet(0.0..0.0))
        "[]"    => intersect(IntervalSet(1.0..1.0), IntervalSet(0.0..0.0))
        "[)"    => intersect(IntervalSet(1.0..1.0), IntervalSet(0.0..0.0))
        "(]"    => intersect(IntervalSet(1.0..1.0), IntervalSet(0.0..0.0))
        "{}"    => intersect(IntervalSet(1.0..1.0), IntervalSet(0.0..0.0))
        "⊘"     => intersect(IntervalSet(1.0..1.0), IntervalSet(0.0..0.0))
        _       => IntervalSet(parse(Interval{float(Int)}, replace(expr, "∞" => "Inf")))
    end
end

"""
    zero_center_interval(val; LB = Closed, RB = Closed) -> IntervalSet
    zero_center_interval(val::Missing; LB = Closed, RB = Closed) = missing
    zero_center_interval(val::IntervalSet; LB = Closed, RB = Closed) = val

Construct a real interval from `-abs(val)` to `abs(val)`.
"""
function zero_center_interval(val; LB = Closed, RB = Closed)
    lb, ub = val > 0 ? (-val, val) : (val, -val)
    if isinf(val)
        LB = Open 
        RB = Open 
    end
    IntervalSet(Interval{float(Int), LB, RB}(lb, ub))
end

zero_center_interval(val::Missing; LB = Closed, RB = Closed) = missing
zero_center_interval(val::IntervalSet; LB = Closed, RB = Closed) = val

# real_interval(center, delta; LB = Closed, RB = Closed) = zero_center_interval(delta; LB, RB) + center

"""
    makecrit_value(crit::Criteria, x) -> Tuple

Create a tuple of criteria based on input value `x`. The relative criterion is multiplied by `x`. 

# Example 
julia> c = crit(10, 0.1) 
Criteria{Int64, Float64}(10, 0.1)

julia> makecrit_value(c, 50)
(10, 5.0)

julia> c = crit(ri"[50,100]", ri"[0.8,1.25]")
Criteria{IntervalSet, IntervalSet}([50.0, 100.0], [0.8, 1.25])

julia> makecrit_value(c, 50)
([50.0, 100.0], [40.0, 62.5])

"""
makecrit_value(c::Criteria, x) = (c.aval, c.rval * x)
makecrit_value(c::Criteria{Missing}, x) = (c.rval * x, )
makecrit_value(c::Criteria{A, Missing}, x) where A = (c.aval, )

"""
    makecrit_delta(crit::Criteria, x) -> Tuple

Create a tuple of criteria based on input value `x`. The relative criterion is multiplied by `x`, and both criteria are offset by `x` as they represent difference.

# Example 
julia> c = crit(10, 0.1) 
Criteria{Int64, Float64}(10, 0.1)

julia> makecrit_delta(c, 50)
([-40.0, 40.0], [-45.0, 45.0])

julia> makecrit_delta(c, 100)
([-90.0, 90.0], [-90.0, 90.0])

"""
makecrit_delta(c::Criteria, x) = (zero_center_interval(c.aval) + x, zero_center_interval(c.rval * x) + x)
makecrit_delta(c::Criteria{Missing}, x) = (zero_center_interval(c.rval * x) + x, )
makecrit_delta(c::Criteria{A, Missing}, x) where A = (zero_center_interval(c.aval) + x, )

_convert(::Type{T}, x::S) where {T, S} = T(x)
_convert(::Type{T}, x::T) where T = x
_convert(::Type{T}, x::Tuple) where T = T.(x)
_convert(::Type{T}, x::NTuple{N, S}) where {T, N, S} = T.(x)
_convert(::Type{T}, x::AbstractVector{S}) where {T, S} = T.(x)
_convert(::Type{T}, x::NTuple{N, T}) where {T, N} = x
_convert(::Type{T}, x::AbstractVector{T}) where T = x
function _vec_macth(regex, p) 
    if eltype(p) <: AbstractString 
        map(p) do x 
            m = match(regex, x)
            isnothing(m) ? 0 : parse(Int, first(m.captures))
        end
    else
        map(p) do x 
            m = match(regex, string(x))
            isnothing(m) ? 0 : parse(Int, first(m.captures))
        end
    end
end

"""
    lastcolnum([output::Type,] ptnames, col; error = true, init = nothing)

Return the last element from `ptnames` such that it starts with `col` and is followed by a number. 
"""
lastcolnum(p, col; error = true, init = nothing) = lastcolnum(Symbol, p, col; error, init)
function lastcolnum(::Type{T}, p, col; error = true, init = nothing) where T
    regex = Regex(string("^", col, "(\\d+\$)"))
    istage = _vec_macth(regex, p)
    stage, icol = findmax(istage)
    stage == 0 && (error ? throw(ArgumentError("No column `$(string(col))...` in table.")) : return init)
    _convert(T, p[icol])
end

"""
    ithcolnum([output::Type,] ptnames, col, i::Int; error = true, init = nothing)

Return element from `ptnames` such that it starts with `col` and is followed by the number `i`. 
"""
ithcolnum(p, col, i::Int; error = true, init = nothing) = ithcolnum(Symbol, p, col, i; error, init)
function ithcolnum(::Type{T}, p, col, i::Int; error = true, init = nothing) where T
    regex = Regex(string("^", col, "(\\d+\$)"))
    istage = _vec_macth(regex, p)
    icol = findfirst(==(i), istage)
    isnothing(icol) && (error ? throw(ArgumentError("No column `$(string(col))$i` in table.")) : return init)
    _convert(T, p[icol])
end

"""
    allcolnum([output::Type,] ptnames, col; error = true, init = [])

Return all elements from `ptnames` such that they start with `col` and is followed by a number. 
"""
allcolnum(p, col; error = true, init = []) = allcolnum(Symbol, p, col; error, init)
function allcolnum(::Type{T}, p, col; error = true, init = []) where T
    regex = Regex(string("^", col, "(\\d+\$)"))
    istage = _vec_macth(regex, p)
    icol = findall(>(0), istage)
    isempty(icol) && (error ? throw(ArgumentError("No column `$(string(col))...` in table.")) : return init)
    id = sortperm(collect(istage[icol]))
    _convert(T, p[icol][id])
end

preabtype(::Val{x}) where x = x
preabtype(::Val{:list}) = :total
preabundance(::Val{x}, abundance, threshold) where x = abundance
function preabundance(::Val{:list}, abundance, threshold) 
    # approximation by [13C] isotope
    m = minimum(makecrit_value(crit(threshold), abundance))
    abundance * (2 * abundance - m) / (abundance - m)
end
dopostnormalize(::Val{x}) where x = false 
dopostnormalize(::Val{:list}) = true
function msmspremaxabundance(vab::Val{x}, abundance, threshold, max_proportion_vec, element_dictionary_vec) where x 
    preab = preabundance(vab, abundance, threshold)
    maxab = preab
    for p in max_proportion_vec
        maxab *= p 
    end
    preab, maxab
end

function msmspremaxabundance(vab::Val{:max}, abundance, threshold, max_proportion_vec, element_dictionary_vec) 
    preab = preabundance(vab, abundance, threshold)
    maxab = preab 
    for p in max_proportion_vec
        preab /= p 
    end
    preab, maxab
end

function msmspremaxabundance(vab::Val{:input}, abundance, threshold, max_proportion_vec, element_dictionary_vec) 
    preab = preabundance(vab, abundance, threshold)
    maxab = preab 
    for p in element_dictionary_vec
        preab /= isotopicabundance(p)
    end
    maxab = preab
    for p in max_proportion_vec
        maxab *= p 
    end
    preab, maxab
end

abtypeop(::Val{:max}) = maximum
abtypeop(::Val{:input}) = first
abtypeop(::Val{:list}) = sum
abtypeop(::Val{:raw}) = nothing
abtypeop(::Val{:total}) = nothing

normab_doc = """
    normalize_abundance(abvector, abundance, abtype, available_abtype = [:max, :input, :list, :raw, :total])
    normalize_abundance!(abvector, abundance, abtype, available_abtype = [:max, :input, :list, :raw, :total])

Normalize raw abundance `abvector` such that the value specified by `abtype` is `abundance`.
"""
@doc normab_doc
normalize_abundance(abvector, abundance, abtype, available_abtype = [:max, :input, :list, :raw, :total]) = 
    _normalize_abundance!(copy(abvector), abundance, abtype, available_abtype; conversion = true)

@doc normab_doc
normalize_abundance!(abvector, abundance, abtype, available_abtype = [:max, :input, :list, :raw, :total]) = 
    _normalize_abundance!(abvector, abundance, abtype, available_abtype)
function _normalize_abundance!(abvector, abundance, abtype, available_abtype = [:max, :input, :list, :raw, :total]; conversion = false)
    isempty(abvector) && return abvector
    abtype in available_abtype || throw(ArgumentError("$abtype is not valid; please select from $available_abtype."))
    fn = abtypeop(Val(abtype))
    isnothing(fn) && return abvector
    abvector = vectorize(abvector)
    T = eltype(abvector)
    S = promote_type(typeof(first(abvector) * abundance / fn(abvector)), T)
    if T == S
        abvector .*= abundance / fn(abvector)
    elseif conversion 
        abvector .* (abundance / fn(abvector))
    else
        throw(ArgumentError("Cannot convert type $S to type $T. Use `normalize_abundance` instead."))
    end
end

tuplize(x::Tuple) = x
tuplize(x) = (x, )
vectorize(x::AbstractDict) = x
vectorize(x::SplitApplyCombine.Dictionaries.AbstractDictionary) = x
vectorize(x::AbstractVector) = x
# vectorize(x::AbstractVector, n::Int) = n == lastindex(x) ? x : n > lastindex(x) ? vcat(x, repeat([last(x)], n - lastindex(x))) : x[begin:begin + n - 1]
vectorize(x) = [x]
vectorize(x, n) = repeat([x], n)