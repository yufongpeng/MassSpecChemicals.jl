"""
    Criteria{A, R}

A type for creating criteria of numeric value.

# Fields
* `aval`: criteria for the absolute value (unnormalized).
* `rval`: criteria for the relative value (normalized by a input)

# Examples
```julia
julia> peak_crit = Criteria(10, 0.2);

julia> qualified_peak(x, basepeak, crit) = all(<=(x), crit(basepeak)); # filter out small peaks based on base peak. 

julia> peak_crit(100) # create criteria a peak value 100
(10, 20)

julia> peak_crit(10) # create criteria a peak value 10
(10, 2)

julia> qualified_peak(10, 100, peak_crit)
false

julia> qualified_peak(20, 100, peak_crit)
true

julia> qualified_peak(10, 10, peak_crit)
true

julia> qualified_peak(2, 10, peak_crit)
false

"""
struct Criteria{A, R}
    aval::A
    rval::R
end

(c::Criteria)(x) = (c.aval, c.rval * x)
(c::Criteria{Missing})(x) = (c.rval * x, )
(c::Criteria{A, Missing})(x) where A = (c.aval, )

abstract type RealIntervals end
struct EmptyInterval <: RealIntervals end
"""
    RealInterval <: RealIntervals

A type representing a real interval. See `@ri_str` for creating a real interval.
"""
struct RealInterval <: RealIntervals
    lowerbound
    upperbound
    leftoperator
    rightoperator
end
"""
    UnionInterval{N} <: RealIntervals

A type representing an union of multiple disjoint real intervals.
"""
struct UnionInterval{N} <: RealIntervals
    intervals::NTuple{N, RealInterval}
end