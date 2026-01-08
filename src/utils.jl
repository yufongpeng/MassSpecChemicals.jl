"""
    acrit(x)

Cnstructing a `Criteria` with only absolute criteria.
"""
acrit(x) = Criteria(x, missing)
"""
    rcrit(x)

Cnstructing a `Criteria` with only relative criteria. 
"""
rcrit(x) = Criteria(missing, x)
"""
    crit(absolute)
    crit(x::Criteria)
    crit(absolute, relative)

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
        _       => IntervalSet(parse(Interval{Float64}, replace(expr, "∞" => "Inf")))
    end
end

"""
    zero_center_interval(val; LB = Closed, RB = Closed)
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
    IntervalSet(Interval{Float64, LB, RB}(lb, ub))
end

zero_center_interval(val::Missing; LB = Closed, RB = Closed) = missing
zero_center_interval(val::IntervalSet; LB = Closed, RB = Closed) = val

# real_interval(center, delta; LB = Closed, RB = Closed) = zero_center_interval(delta; LB, RB) + center

"""
    makecrit_value(crit, x)

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
    makecrit_delta(crit, x)

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

function makerttol(fwhm_ion, fwhm_lib, rrt, lib)
    fwhm = isnothing(fwhm_lib) ? vectorize(fwhm_ion) : [(fwhm_ion + x) / 2 for x in vectorize(fwhm_lib)]
    if length(fwhm) == 1
        [union(makecrit_delta(crit(first(fwhm)), rrt)...) for _ in eachindex(lib)]
    elseif length(fwhm) != length(lib)
        [union(makecrit_delta(crit(mean(fwhm)), rrt)...) for _ in eachindex(lib)]
    else
        [union(makecrit_delta(crit(x), rrt)...) for x in fwhm]
    end
end
function collecfwhm(fwhm, exp)
    v = vectorize(fwhm)
    if length(v) == 1
        repeat(v, length(exp))
    elseif length(v) != length(exp)
        repeat(mean(v), length(exp))
    else
        v
    end
end

tuplize(x::Tuple) = x
tuplize(x) = (x, )
vectorize(x::AbstractVector) = x
vectorize(x::AbstractVector, n::Int) = x[begin:begin + n - 1]
vectorize(x) = [x]
vectorize(x, n) = repeat([x], n)