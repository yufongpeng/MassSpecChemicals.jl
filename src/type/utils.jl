"""
    Criteria{A, R}

A type for creating criteria of numeric value. Criterion may be a cutoff value, allowed tolerance relative to an input value, or acceptable range (interval(s)).

# Fields
* `aval`: criterion for the absolute value (unnormalized).
* `rval`: criterion for the relative value (normalized by an input)

# Examples
```julia
julia> peak_crit = Criteria(10, 0.2);

julia> qualified_peak(x, basepeak, crit) = all(<=(x), makecrit_value(crit, basepeak)); # filter out small peaks based on base peak. 

julia> makecrit_value(peak_crit, 100) # create criteria for peak value 100
(10, 20)

julia> makecrit_value(peak_crit, 10) # create criteria for peak value 10
(10, 2)

julia> qualified_peak(10, 100, peak_crit)
false

julia> qualified_peak(20, 100, peak_crit)
true

julia> qualified_peak(10, 10, peak_crit)
true

julia> qualified_peak(2, 10, peak_crit)
false

julia> qc_peak(x, refpeak, crit) = all(in(x), makecrit_delta(crit, refpeak)); # filter out deviated peaks.

julia> makecrit_delta(peak_crit, 100) # create criteria for reference peak value 100
(Intervals.Interval{Float64, Intervals.Closed, Intervals.Closed}(90.0, 110.0), Intervals.Interval{Float64, Intervals.Closed, Intervals.Closed}(80.0, 120.0))

julia> qc_peak(121, 100, peak_crit)
false

"""
struct Criteria{A, R}
    aval::A
    rval::R
end

"""
    value_error(x, y)

Error function of difference; `y - x`.
"""
value_error(x, y) = y - x

"""
    relative_error(x, y)

Error function of relative difference (relative to true value x); `(y - x) / x`.
"""
relative_error(x, y) = (y - x) / x

"""
    relative_error_mean(x, y)

Error function of relative difference (relative to mean); `(y - x) / ((x + y) / 2)`.
"""
relative_error_mean(x, y) = 2 * (y - x) / (x + y) 

"""
    percentage_error(x, y)

Error function of relative difference in percentage (relative to true value x); `(y - x) / x * 100`.
"""
percentage_error(x, y) = (y - x) / x * 100

"""
    percentage_error_mean(x, y)

Error function of relative difference in percentage (relative to mean); `(y - x) / ((x + y) / 2) * 100`.
"""
percentage_error_mean(x, y) = 2 * (y - x) / (x + y) * 100

"""
    ppm_error(x, y)

Error function of relative difference in ppm (relative to true value x); `(y - x) / x * 1e6`.
"""
ppm_error(x, y) = (y - x) / x * 1e6

"""
    ppm_error_mean(x, y)

Error function of relative difference in ppm (relative to mean); `(y - x) / ((x + y) / 2) * 1e6`.
"""
ppm_error_mean(x, y) = 2 * (y - x) / (x + y) * 1e6