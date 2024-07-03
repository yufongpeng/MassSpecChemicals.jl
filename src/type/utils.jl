struct Criteria{A, R}
    aval::A
    rval::R
end

abstract type RealIntervals end
struct EmptyInterval <: RealIntervals end
"""
    RealInterval <: RealIntervals

A type representing a real interval.
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