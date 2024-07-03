*(x::Criteria, y::Number) = Criteria(x.aval * y, x.rval * y)
/(x::Criteria, y::Number) = Criteria(x.aval / y, x.rval / y)

*(x::RealInterval, y::Number) = RealInterval(x.lowerbound * y, x.upperbound * y, x.leftoperator, x.rightoperator)
/(x::RealInterval, y::Number) = RealInterval(x.lowerbound / y, x.upperbound / y, x.leftoperator, x.rightoperator)

function isless(x::RealInterval, y::RealInterval)
    x.lowerbound < y.lowerbound && return true
    x.lowerbound > y.lowerbound && return false
    xor((x.leftoperator)(x.lowerbound, y.lowerbound), (y.leftoperator)(x.lowerbound, y.lowerbound)) && return (x.leftoperator)(x.lowerbound, y.lowerbound)
    x.upperbound < y.upperbound && return true
    x.upperbound > y.upperbound && return false
    xor((x.rightoperator)(x.upperbound, y.upperbound), (y.rightoperator)(x.upperbound, y.upperbound)) && return (y.rightoperator)(x.upperbound, y.upperbound)
    false
end

isequal(x::RealInterval, y::RealInterval) = 
    x.lowerbound == y.lowerbound && x.upperbound == y.upperbound && !xor((x.rightoperator)(x.upperbound, y.upperbound), (y.rightoperator)(x.upperbound, y.upperbound))

in(x::Number, ri::EmptyInterval) = false
in(x::Number, ri::RealInterval) = between(x; low = ri.lowerbound, up = ri.upperbound, lop = ri.leftoperator, rop = ri.rightoperator)
in(x::Number, ui::UnionInterval) = any(ri -> between(x; low = ri.lowerbound, up = ri.upperbound, lop = ri.leftoperator, rop = ri.rightoperator), ui.intervals)

function union(ris::Vararg{<: RealIntervals, N}) where N
    vri = RealInterval[]
    for ri in ris
        @match ri begin
            ::UnionInterval => union!(vri, ri.intervals)
            ::RealInterval  => push!(vri, ri)
            ::EmptyInterval => nothing
        end
    end
    unique!(vri)
    sort!(vri)
    nvri = RealInterval[]
    i = 1
    push!(nvri, vri[1])
    while i < length(vri)
        i += 1
        if nvri[end].upperbound in vri[i]
            nvri[end] = RealInterval(nvri[end].lowerbound, vri[i].upperbound, nvri[end].leftoperator, vri[i].rightoperator)
        elseif vri[i].lowerbound in nvri[end]
            if nvri[end].upperbound == vri[i].lowerbound
                nvri[end] = RealInterval(nvri[end].lowerbound, vri[i].upperbound, nvri[end].leftoperator, vri[i].rightoperator)
            end
        else
            push!(nvri, vri[i])
        end
    end
    length(nvri) == 1 ? nvri[1] : UnionInterval(Tuple(nvri))
end

intersect(ris::Vararg{<: Union{RealInterval, UnionInterval}, N}) where N = reduce(intersect, ris)
intersect(ri1::EmptyInterval, ri2::RealIntervals) = EmptyInterval()
intersect(ri1::RealIntervals, ri2::EmptyInterval) = EmptyInterval()
intersect(ri1::EmptyInterval, ri2::EmptyInterval) = EmptyInterval()

function intersect(ri1::RealInterval, ri2::RealInterval)
    ri1 === ri2 && return ri1
    if ri1 > ri2
        ri1, ri2 = ri2, ri1
    end
    if !xor((ri2.rightoperator)(ri1.upperbound, ri2.upperbound), (ri1.rightoperator)(ri2.upperbound, ri1.upperbound))
        return ri2
    elseif ri1.upperbound == ri2.upperbound
        return (ri2.rightoperator)(ri1.upperbound, ri2.upperbound) ? RealInterval(ri2.lowerbound, ri1.upperbound, ri2.leftoperator, ri1.rightoperator) : ri2
    elseif (ri1.rightoperator)(ri2.upperbound, ri1.upperbound)
        return ri2
    elseif (ri2.leftoperator)(ri2.lowerbound, ri1.upperbound) && (ri1.rightoperator)(ri2.lowerbound, ri1.upperbound)
        return RealInterval(ri2.lowerbound, ri1.upperbound, ri2.leftoperator, ri1.rightoperator)
    else 
        return EmptyInterval()
    end
end

intersect(ri1::UnionInterval, ri2::RealInterval) = union((intersect(ri, ri2) for ri in ri1)...)
intersect(ri1::RealInterval, ri2::UnionInterval) = union((intersect(ri1, ri) for ri in ri2)...)
intersect(ri1::UnionInterval, ri2::UnionInterval) = union((intersect(ria, rib) for ria in ri1, rib in ri2)...)
