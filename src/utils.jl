(c::Criteria)(x) = (c.aval, c.rval * x)
(c::Criteria{Missing})(x) = (c.rval * x, )
(c::Criteria{A, Missing})(x) where A = (c.aval, )

acrit(x) = Criteria(x, missing)
rcrit(x) = Criteria(missing, x)
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
    @ri_str -> RealInterval

Create a `RealInterval` by mathematical real interval notation.

# Examples
```julia
julia> f = ri"[4, 7)"
(::MSChemicals.RealInterval) (generic function with 1 method)

julia> f(4)
true

julia> f(7)
false
```
"""
macro ri_str(expr)
    return _real_interval(expr)
end
function _real_interval(expr::AbstractString)
    lc, lv, rv, rc = match(r" *([\(\[]) *([+-]*[\d∞Inf]*\.*\d*) *, *([+-]*[\d∞Inf]*\.*\d*) *([\)\]]) *", expr)
    lop = @match lc begin
        "[" => <=
        "(" => <
    end
    rop = @match rc begin
        "]" => <=
        ")" => <
    end
    lv = @match lv begin
        "+∞"    => Inf
        "+Inf"  => Inf
        "∞"     => Inf
        "Inf"   => Inf
        "-∞"    => -Inf
        "-Inf"  => -Inf
        if occursin(".", lv) end   => parse(Float64, lv)
        _       => parse(Int, lv)
    end
    rv = @match rv begin
        "+∞"    => Inf
        "+Inf"  => Inf
        "∞"     => Inf
        "Inf"   => Inf
        "-∞"    => -Inf
        "-Inf"  => -Inf
        if occursin(".", rv) end   => parse(Float64, rv)
        _       => parse(Int, rv)
    end
    _real_interval(lv, rv, lop, rop)
end

function _real_interval(lb, ub, lop = <=, rop = <=)
    lop(lb, ub) || return EmptyInterval()
    rop(lb, ub) || return EmptyInterval()
    RealInterval(lb, ub, lop, rop)
end

function real_interval(val, lop = <=, rop = <=)
    lb, ub = val > 0 ? (-val, val) : (val, -val)
    lop(lb, ub) || return EmptyInterval()
    rop(lb, ub) || return EmptyInterval()
    RealInterval(lb, ub, lop, rop)
end

real_interval(val::Missing, lop = <=, rop = <=) = missing
real_interval(val::RealIntervals, lop = <=, rop = <=) = val

function real_interval(ct::Criteria, lop = <=, rop = <=)
    Criteria(real_interval(ct.aval, lop, rop), real_interval(ct.rval, lop, rop))
end

tuplize(x::Tuple) = x
tuplize(x) = (x, )
vectorize(x::AbstractVector) = x
vectorize(x::AbstractVector, n::Int) = x[begin:begin + n]
vectorize(x) = [x]
vectorize(x, n) = repeat([x], n)