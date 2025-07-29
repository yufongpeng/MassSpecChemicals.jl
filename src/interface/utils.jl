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