"""
    MSScan([msanalyzer = TOF(),] mztable; min_bin_fwhm = 50) -> Spectrum
    MSScan([msanalyzer = TOF(),] spectrum; min_bin_fwhm = 50) -> Spectrum

Perform MS scan. Theoretical m/z signals are convoluted with `msanalyzer.window`. 

* `msanalyzer::AbstractMSAnalyzer`: a MS analyzer to perform MS Scan. See documentation of specific analyzer for detailed settings.
* `mztable::Table`: a table containing columns
    * `Abundance1`, `Abundance2`, ..., `Abundancen`. The last column will be utilized.
    * `MZ1`, `MZ2`, ..., `MZn`. The last column will be utilized.
* `spectrum::Spectrum`: `spectrum.table` is utilized as `mztable`.
* `min_bin_fwhm`: minimal number of bins within fwhm. 
"""
function MSScan(msanalyzer::AbstractMSAnalyzer, mztable::Table; min_bin_fwhm = 50)
    digits = hasproperty(msanalyzer, :digits) ? msanalyzer.digits : nothing
    stepsize = hasproperty(msanalyzer, :stepsize) ? msanalyzer.stepsize : nothing
    if !isnothing(stepsize) && !isnothing(digits) 
        unit = round(10.0 ^ (-digits); digits)
        p = round(Int, stepsize / unit)
        isapprox(stepsize - p * unit, 0) || throw(ArgumentError("`stepsize` must be mutiples of `10 ^ (-digits)`."))
    end
    icol = findlast(x -> startswith(x, "MZ"), string.(propertynames(mztable)))
    isnothing(icol) && throw(ArgumentError("No column `MZ1`, ..., `MZn` in mztable."))
    colmz = propertynames(mztable)[icol]
    icol = findlast(x -> startswith(x, "Abundance"), string.(propertynames(mztable)))
    isnothing(icol) && throw(ArgumentError("No column `Abundance1`, ..., `Abundacnen` in mztable."))
    colab = propertynames(mztable)[icol]
    id = sortperm(getproperty(mztable, colmz))
    if !isnothing(msanalyzer.mz)
        lower_mz, upper_mz = msanalyzer.mz
        mz_vector = getproperty(mztable, colmz)
        filter!(x -> lower_mz < mz_vector[x] < upper_mz, id)
    end
    mztable = mztable[id]
    mz_vector = getproperty(mztable, colmz)
    ab_vector = getproperty(mztable, colab)
    fwhm = fwhm_mz(msanalyzer, first(mz_vector))
    binsize = min(makecrit_value(crit(msanalyzer.accuracy), first(mz_vector))..., fwhm / min_bin_fwhm)
    nbin_multiplier = 1
    if !isnothing(stepsize)
        nbin_multiplier = ceil(Int, stepsize / binsize)
        binsize = stepsize / nbin_multiplier
    end
    kernels = discrete_window(msanalyzer.window, mz_vector, msanalyzer, binsize, nbin_multiplier, minimum(ab_vector) / maximum(ab_vector) / 10)
    # ihwhm, k = convolution_window(window, hwhm, binsize, nbin_multiplier, minimum(mztable.Abundance) / 10; normalize = true)
    initmass = isnothing(digits) ? first(mz_vector) : round(first(mz_vector); digits)
    binmass, ibins = binnify(mz_vector, binsize, initmass)
    spectrum_offset = bin_offset(initmass, first(binmass), binsize, nbin_multiplier)
    convolution_vector = [a .* k for (a, k) in zip(ab_vector, kernels)]
    ibins .+= Int(length(first(kernels))) ÷ 2 + 1
    spectrum = zeros(eltype(ab_vector), last(ibins) + Int(length(last(kernels))) ÷ 2)
    for (c, i) in zip(convolution_vector, ibins)
        convolution_offset = Int(length(c) ÷ 2)
        spectrum[i - convolution_offset:i + convolution_offset] .+= c 
    end
    # spectrum = spectrum_offset > 1 ? spectrum[spectrum_offset:end] : spectrum
    spectrum = spectrum_offset > 0 ? vcat(zeros(spectrum_offset), spectrum) : spectrum
    # ibins .-= spectrum_offset - 1
    ibins .+= spectrum_offset
    mztable = Table(mztable; 
        Bin_MZ = binmass, 
        Convolution = convolution_vector, 
        Bin_id = ibins
    )
    initial_mass = first(binmass) - (first(ibins) - 1) * binsize
    # pt = peak_table(mztable, spectrum, initial_mass, binsize)
    Spectrum(spectrum, initial_mass, binsize, nbin_multiplier, mztable)
end
MSScan(mztable; min_bin_fwhm = 50) = MSScan(TOF(), mztable; min_bin_fwhm)
MSScan(msanalyzer::AbstractMSAnalyzer, spec::Spectrum; min_bin_fwhm = 50) = MSScan(msanalyzer, spec.table; min_bin_fwhm)

"""
    AllIons([mz_range = nothing,] mztable) -> Table
    AllIons([mz_range = nothing,] spectrum) -> Table

Allow all Ions within m/z range entering the next MS stage. 

* `mz_range::Union{Nothing, Tuple}`: nothing (indicating all ions) or a tuple of m/z lower bound an d upper bound.
* `mztable::Table`: a table containing columns
    * `MZ1`, `MZ2`, ..., `MZn`. The last column will be utilized.
* `spectrum::Spectrum`: `spectrum.table` is utilized as `mztable`.
"""
function AllIons(mz_range, mztable::Table)
    icol = findlast(x -> startswith(x, "MZ"), string.(propertynames(mztable)))
    isnothing(icol) && throw(ArgumentError("No column `MZ1`, ..., `MZn` in mztable."))
    if !isnothing(mz_range)
        lower_mz, upper_mz = mz_range 
        id = findall(x -> lower_mz < x < upper_mz, getproperty(mztable, propertynames(mztable)[icol]))
        mztable = mztable[id]
    end
    mztable
end
AllIons(mztable) = AllIons(nothing, mztable)
AllIons(mz_range, spec::Spectrum) = AllIons(mz_range, spec.table)

"""
    TargetIon(msanalyzer, mztable; threshold = rcrit(1e-4)) -> Table
    TargetIon(msanalyzer, spectrum; threshold = rcrit(1e-4)) -> Table

Isolating target ion(s) with specific m/z values and resolution to enter the next MS stage. 

* `msanalyzer::AbstractMSAnalyzer`: a MS analyzer to perform MS filtering. See documentation of specific analyzer for detailed settings.
* `mztable::Table`: a table containing columns    
    * `ID`: ID tuples. Each elements represents ID number of ions of each MS stage. 
    * `Abundance1`, `Abundance2`, ..., `Abundancen`. The last column will be utilized.
    * `MZ1`, `MZ2`, ..., `MZn`. The last column will be utilized.
* `spectrum::Spectrum`: `spectrum.table` is utilized as `mztable`.
* `threshold` can be a number or criteria (absolute and/or relative to maximum), representing the lower limit of abundance. 
"""
function TargetIon(msanalyzer::AbstractMSAnalyzer{W, M}, mztable::Table; threshold = rcrit(1e-4)) where {W, M <: Real}
    icol = findlast(x -> startswith(x, "MZ"), string.(propertynames(mztable)))
    isnothing(icol) && throw(ArgumentError("No column `MZ1`, ..., `MZn` in mztable."))
    colmz = propertynames(mztable)[icol] 
    icol = findlast(x -> startswith(x, "Abundance"), string.(propertynames(mztable)))
    isnothing(icol) && throw(ArgumentError("No column `Abundance1`, ..., `Abundacnen` in mztable."))
    colab = propertynames(mztable)[icol]
    params = window_parameter(msanalyzer, msanalyzer.mz)
    ab = map(mztable) do r
        getproperty(r, colab) * msanalyzer.window(getproperty(r, colmz), params...)
    end
    tab = minimum(makecrit_value(crit(threshold), maximum(ab)))
    id = findall(>(tab), ab)
    Table(mztable; [colab => ab]...)[id]
end
TargetIon(msanalyzer::AbstractMSAnalyzer{W, M}, spec::Spectrum; threshold = rcrit(1e-4)) where {W, M <: Real} = TargetIon(msanalyzer, spec.table; threshold)

"""
    Fragmentation(product_table, precursor_table; threshold) -> Table
    Fragmentation(product_table, spectrum; threshold) -> Table

Fragmentation of `precursor_table.Chemical` or `spectrum.table.Chemical` into `product_table.Product`.

* `producttable::Table`: a table containing columns
    * `ID`: ID tuples map product information to precursors in `precursor_table.Chemical`.
    * `Product`: products of each precursor. 
    * `Proportion`: proportion of fragmentation of each product relative to precursor signal. This column is optional; the default is that each product share precursor signals equally. 
* `precursortable::Table`: a table containing columns
    * `ID`: ID tuples. Each element represents ID number of ions of each MS stage. 
    * `Chemical`: chemical objects. 
    * `Parent`: parent chemical objects. This column is required if `Chemical`s are formulas.
    * `Abundance1`, `Abundance2`, ..., `Abundancen`. The last column will be utilized.
    * `MZ1`, `MZ2`, ..., `MZn`. The last column will be utilized.
* `spectrum::Spectrum`: `spectrum.table` is utilized as `precursortable`.
* `threshold` can be a number or criteria, representing the lower limit of abundance (absolute and/or relative to maximal value of each spectrum). 
"""
function Fragmentation(producttable::Table, mztable::Table; threshold = rcrit(1e-4))
    # group -> precursor_table, elements_precursor
    if !in(:ID, propertynames(producttable)) && in(:Chemical, propertynames(producttable))
        producttable = match_chemical(mztable, producttable; colexp = :Chemical, collib = :Chemical)
    else
        throw(ArgumentError("No column `ID` or `Chemical` in product_table."))
    end
    :ID in propertynames(mztable) || throw(ArgumentError("No column `ID` in precursor_table."))
    :Chemical in propertynames(mztable) || throw(ArgumentError("No column `Chemical` in precursor_table."))
    :Abundance1 in propertynames(mztable) || throw(ArgumentError("No column `Abundance1`, ..., `Abundacnen` in precursor_table."))
    :MZ1 in propertynames(mztable) || throw(ArgumentError("No column `MZ1`, ..., `MZn` in precursor_table."))
    :Product in propertynames(producttable) || throw(ArgumentError("No column `Product` in product_table."))
    allequal(msstage, mztable.Chemical) || throw(ArgumentError("Chemicals have to be in the same MS stage."))
    all(x -> all(y -> msstage(y) < 2, x), producttable.Product) || throw(ArgumentError("Products should not MS/MS pairs."))
    if !in(:Proportion, propertynames(producttable))
        producttable = Table(producttable; Proportion = [nothing for _ in eachindex(producttable)])
    end
    # MS1 threshold
    # threshold = acrit(minimum(makecrit_value(crit(threshold), maximum(mztable.Abundance1))))
    id = producttable.ID
    mztable = filter(x -> x.ID in id, mztable)
    gmztable = group(getproperty(:ID), mztable)
    vcat(map(gmztable) do precursor_table 
        pid = findfirst(==(first(precursor_table.ID)), id)
        TandemIsotopologues(chemicalparent(first(precursor_table.Chemical)); 
            threshold,
            precursor_table,
            product = producttable.Product[pid],
            proportion = producttable.Proportion[pid]
            )
    end...)
end
Fragmentation(producttable::Table, spec::Spectrum; threshold = rcrit(1e-4)) = Fragmentation(producttable, spec.table; threshold)

"""
    peak_table(spectrum::Spectrum; alg = LocalMaxima(), abundance = 1, abtype = :max, threshold = crit(abundance * 1e-4, 1e-4))

Extract peaks from a spectrum.

* `abundance` sets the abundance of the peak specified by `abtype`. 
* `abtype`
    * `:max`: the largest peak.
    * `:list`: sum of listed peaks.
    * `:raw`: no abundance normalization.
* `threshold` can be a number or criteria (absolute and/or relative to `abundance`), representing the lower limit of abundance. 
"""
peak_table(spectrum::Spectrum; alg = LocalMaxima(), abundance = 1, abtype = :max, threshold = rcrit(1e-4)) = peak_table(spectrum.table, spectrum.spectrum, spectrum.initial_mass, spectrum.binsize, spectrum.stepsize; alg, abundance, abtype, threshold)
function peak_table(mztable::Table, spectrum, initial_mass, binsize, nbin_multiplier; alg = LocalMaxima(), abundance = 1, abtype = :max, threshold = rcrit(1e-4))
    maxbin = [find_nearest_peak(alg, spectrum, r.Bin_id, r.Convolution) for r in mztable] 
    if nbin_multiplier > 1
        maxbin = map(maxbin) do mbin 
            isnothing(mbin) && return nothing
            r = rem(mbin - 1, nbin_multiplier)
            nbin_multiplier - r < r ? mbin + nbin_multiplier - r : mbin - r 
        end
    end
    id = findall(!isnothing, maxbin)
    mztable = Table(mztable[id]; Max_bin = maxbin[id])
    gmztable = group(getproperty(:Max_bin), mztable)
    c = length(first(mztable.Convolution)) ÷ 2 + 1
    tuples = map(pairs(gmztable)) do (ibin, smztable) 
        cab = map(smztable) do r 
            delta = ibin - r.Bin_id
            r.Convolution[c + delta]
        end
        id = sortperm(cab; rev = true)
        (; Chemical = Isobars(getproperty(smztable, :Chemical)[id], cab[id]), MZ = initial_mass + (ibin - 1) * binsize, Abundance = spectrum[ibin])
    end
    table = Table([t for t in tuples])
    if abtype == :max
        table.Abundance .*= abundance / maximum(table.Abundance) 
    elseif abtype == :list 
        table.Abundance .*= abundance / sum(table.Abundance) 
    end
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(table.Abundance)))
    table[findall(>=(abundance_cutoff), table.Abundance)]
end

"""
    plot_spectrum([mz_range = nothing,] spectrum::Spectrum; deconvolution = false, threshold = rcrit(1e-4), kwargs...)
    plot_spectrum([mz_range = nothing,] mztable::Table; threshold = rcrit(1e-4), kwargs...)

Plot a spectrum.

* `mz_range::Union{Nothing, Tuple}`: nothing (indicating entire mz range) or a tuple of m/z lower bound an d upper bound.
* `deconvolution::Bool`: plot deconvoluted or convoluted spectrum.
* `threshold` can be a number or criteria (absolute and/or relative to maximum), data lower than this value will be set to zero. 

Other keyword arguments can controls the settings of plot.
"""
plot_spectrum(mz_range, spectrum::Spectrum; deconvolution = false, threshold = rcrit(1e-4), kwargs...) = 
    _plot_spectrum(mz_range, spectrum; deconvolution, threshold, kwargs...)
plot_spectrum(mz_range, mztable::Table; threshold = rcrit(1e-4), kwargs...) = 
    _plot_spectrum(mz_range, mztable; threshold, kwargs...)
plot_spectrum(x; kwargs...) = plot_spectrum(nothing, x; kwargs...)

"""
    plot_spectrum!([mz_range = nothing,] spectrum::Spectrum; deconvolution = false, threshold = rcrit(1e-4), kwargs...)
    plot_spectrum!([mz_range = nothing,] mztable::Table; threshold = rcrit(1e-4), kwargs...)

Plot a spectrum to an existing figure.

* `mz_range::Union{Nothing, Tuple}`: nothing (indicating entire mz range) or a tuple of m/z lower bound an d upper bound.
* `deconvolution::Bool`: plot deconvoluted or convoluted spectrum.
* `threshold` can be a number or criteria (absolute and/or relative to maximum), data lower than this value will be set to zero. 

Other keyword arguments can controls the settings of plot.
"""
plot_spectrum!(mz_range, spectrum::Spectrum; deconvolution = false, threshold = rcrit(1e-4), kwargs...) = 
    _plot_spectrum(mz_range, spectrum; fn = plot!, deconvolution, threshold, kwargs...)
plot_spectrum!(mz_range, mztable::Table; threshold = rcrit(1e-4), kwargs...) = 
    _plot_spectrum(mz_range, mztable; fn = plot!, threshold, kwargs...)
plot_spectrum!(x; kwargs...) = plot_spectrum!(nothing, x; kwargs...)

function _plot_spectrum(mz_range, spectrum::Spectrum; fn = plot, deconvolution = false, threshold = rcrit(1e-4), kwargs...)
    deconvolution && return _plot_spectrum(mz_range, spectrum.table; fn, threshold, kwargs...)
    np =  first(Plots.default(:size))
    if isnothing(mz_range) 
        range = length(spectrum.spectrum) * spectrum.binsize
        mz_lower = max(spectrum.initial_mass - 0.1 * range, 0)
        mz_upper = spectrum.initial_mass + 1.1 * range
    else
        mz_lower, mz_upper = mz_range 
    end
    mzid_s = Int((mz_lower - spectrum.initial_mass) ÷ (spectrum.binsize * spectrum.stepsize)) * spectrum.stepsize
    mzid_e = ceil(Int, (mz_upper - spectrum.initial_mass) / (spectrum.binsize * spectrum.stepsize)) * spectrum.stepsize
    fid = max(firstindex(spectrum.spectrum), mzid_s + 1)
    lid = min(lastindex(spectrum.spectrum), mzid_e + 1)
    delta, r = divrem(mzid_e - mzid_s, np)
    fid >= lid && throw(ArgumentError("No spectrum data; please adjust m/z range."))
    ab_cutoff = minimum(makecrit_value(crit(threshold), maximum(@view spectrum.spectrum[fid:lid])))
    delta = max(delta, 1)
    mzid = delta > 1 ? (mzid_s:delta:(mzid_e - r + delta)) : mzid_s:mzid_e
    spec = map(mzid) do i 
        if i + delta < fid || i + 1 > lid
            s = 0 
        elseif i + 1 < fid
            s = maximum(@view spectrum.spectrum[fid:i + delta])
        elseif i + delta > lid
            s = maximum(@view spectrum.spectrum[i + 1:lid])
        else
            s = maximum(@view spectrum.spectrum[i + 1:i + delta])
        end
        s < ab_cutoff ? 0 : s
    end
    kwargs = Dict(kwargs...)
    spec_kwargs!(kwargs)
    ms = @. spectrum.initial_mass + (mzid * spectrum.binsize)
    fn(ms, spec; kwargs...)
end

function _plot_spectrum(mz_range, mztable::Table; fn = plot, threshold = rcrit(1e-4), kwargs...)
    icol = findlast(x -> startswith(x, "MZ"), string.(propertynames(mztable)))
    isnothing(icol) && throw(ArgumentError("No MZ data to process!"))
    colmz = propertynames(mztable)[icol]
    mz = getproperty(mztable, colmz)
    icol = findlast(x -> startswith(x, "Abundance"), string.(propertynames(mztable)))
    isnothing(icol) && throw(ArgumentError("No Abundance data to process!"))
    colab = propertynames(mztable)[icol]
    if isnothing(mz_range) 
        min_mz, max_mz = extrema(mz)
        range = (max_mz - min_mz) * 0.1
        mz_lower = max(min_mz - range, 0)
        mz_upper = max_mz + range
    else
        mz_lower, mz_upper = mz_range 
    end
    kwargs = Dict(kwargs...)
    spec_kwargs!(kwargs)
    id = findall(x -> mz_lower <= x <= mz_upper, mz)
    gmztable = group(getproperty(colmz), mztable[id])
    ab = map(x -> sum(getproperty(x, colab)), gmztable)
    ab_cutoff = minimum(makecrit_value(crit(threshold), maximum(ab)))
    id = findall(>=(ab_cutoff), ab)
    ab = [ab[i] for i in id]
    fn(vcat(mz_lower, repeat(collect(id); inner = 3), mz_upper), vcat(0, [[0, a, 0] for a in ab]..., 0); kwargs...)
end

function spec_kwargs!(kwargs)
    get!(kwargs, :xlabel, "m/z")
    get!(kwargs, :ylabel, "Abundance")
    get!(kwargs, :title, "Spectrum")
    get!(kwargs, :label, nothing)
    kwargs
end

"""
    plot_resolving_power([mz_range = nothing,] ms::AbstractMSAnalyzer; n = 1000, kwargs...)

Plot the function of m/z to resolving_power.

* `mz_range::Union{Nothing, Tuple}`: nothing (indicating using `ms.mz`) or a tuple of m/z lower bound an d upper bound.
* `n::Int`: number of data points.

Other keyword arguments can controls the settings of plot.
"""
plot_resolving_power(mz_range, ms::AbstractMSAnalyzer; n = 1000, kwargs...) = 
    _plot_resolving_power(mz_range, ms; n, kwargs...)
plot_resolving_power(ms::AbstractMSAnalyzer; n = 1000, kwargs...) = plot_resolving_power(nothing, ms; n, kwargs...)

"""
    plot_resolving_power!([mz_range = nothing,] ms::AbstractMSAnalyzer; n = 1000, kwargs...)

Plot the function of m/z to resolving_power to an existing figure.

* `mz_range::Union{Nothing, Tuple}`: nothing (indicating using `ms.mz`) or a tuple of m/z lower bound an d upper bound.
* `n::Int`: number of data points.

Other keyword arguments can controls the settings of plot.
"""
plot_resolving_power!(mz_range, ms::AbstractMSAnalyzer; n = 1000, kwargs...) = 
    _plot_resolving_power(mz_range, ms; fn = plot!, n, kwargs...)
plot_resolving_power!(ms::AbstractMSAnalyzer; n = 1000, kwargs...) = plot_resolving_power!(nothing, ms; n, kwargs...)

function _plot_resolving_power(mz_range, ms::AbstractMSAnalyzer; fn = plot, n = 1000, kwargs...)
    if isnothing(mz_range)
        mz_range = ms.mz
    end
    if mz_range isa Tuple{<: Real, <: Real}
        mz_lower, mz_upper = mz_range 
    else
        throw(ArgumentError("`mz_range` must be a tuple, i.e. (lowerbound, upperbound)."))
    end
    x = mz_lower:((mz_upper - mz_lower) / n):mz_upper
    y = resolving_power.(Ref(ms), x)
    kwargs = Dict(kwargs...)
    respow_kwargs!(kwargs)
    fn(x, y; kwargs...)
end

function respow_kwargs!(kwargs)
    get!(kwargs, :xlabel, "m/z")
    get!(kwargs, :ylabel, "Resolution (m/Δm)")
    get!(kwargs, :title, "Resolving Power")
    get!(kwargs, :label, nothing)
    kwargs
end

"""
    plot_window(window::AbstractWindow; fwhm = 0.7, binsize = 0.01, nbin_multiplier = 1, height = 0.01, kwargs...)

Plot the window function.

* `fwhm::Real`: m/z fwhm.
* `binsize::Real`: binsize of m/z value. 
* `nbin_multiplier`: multiplier of actual plotted data bin relative to `binsize`.
* `height::Real`: the lower limit of transmission. 

Other keyword arguments can controls the settings of plot.
"""
plot_window(window::AbstractWindow; fwhm = 0.7, binsize = 0.01, stepsize = 1, height = 0.01, kwargs...) = 
    _plot_window(window; fwhm, binsize, stepsize, height)

"""
    plot_window!(window::AbstractWindow; fwhm = 0.7, binsize = 0.01, nbin_multiplier = 1, height = 0.01, kwargs...)

Plot the window function to an existing figure.

* `fwhm::Real`: m/z fwhm.
* `binsize::Real`: binsize of m/z value. 
* `nbin_multiplier`: multiplier of actual plotted data bin relative to `binsize`.
* `height::Real`: the lower limit of transmission. 

Other keyword arguments can controls the settings of plot.
"""
plot_window!(window::AbstractWindow; fwhm = 0.7, binsize = 0.01, stepsize = 1, height = 0.01, kwargs...) = 
    _plot_window(window; fn = plot!, fwhm, binsize, stepsize, height)

function _plot_window(window::AbstractWindow; fn = plot, fwhm = 0.7, binsize = 0.01, stepsize = 1, height = 0.01, kwargs...)
    y = discrete_window(window, fwhm, binsize, stepsize, height)
    x = (eachindex(y) .- (length(y) ÷ 2 + 1)) .* binsize
    kwargs = Dict(kwargs...)
    window_kwargs!(kwargs)
    fn(x, y; kwargs...)
end

function window_kwargs!(kwargs)
    get!(kwargs, :xlabel, "m/z")
    get!(kwargs, :ylabel, "Transmission")
    get!(kwargs, :title, "Window")
    get!(kwargs, :label, nothing)
    kwargs
end

# function convolution_window(window::AbstractWindow, hwhm::T, binsize::S, nbin_multiplier::Int, height; normalize = false) where {T, S}
#     R = promote_type(T, S)
#     hwhm = convert(R, hwhm)
#     binsize = convert(R, binsize)
#     k = discrete_window(window, hwhm, binsize, nbin_multiplier, height)
#     h = first(k) / 2
#     i = findfirst(<(h), k)
#     ihwhm = i - (h - k[i]) / (k[i - 1] - k[i])
#     k = normalize ? k ./ first(k) : k
#     ceil(Int, ihwhm), vcat(reverse(k[begin + 1:end]), k)
# end


function bin_offset(outmass, binmass, binsize, nbin_multiplier)
    x = round(Int, (outmass - binmass) / binsize)
    # x == 0 ? 1 : x > 0 ? x + 1 : nbin_multiplier + x + 1
    x == 0 ? 0 : x > 0 ? nbin_multiplier - x : - x
end

function binnify(mass, binsize, init = first(mass))
    binmass = similar(mass)
    ibins = zeros(Int, length(binmass))
    ibin = 0
    for (i, m) in enumerate(mass)
        while init < m 
            ibin += 1
            init += binsize
        end
        if init - m > binsize / 2 
            binmass[i] = init - binsize 
            ibins[i] = ibin - 1
        else
            binmass[i] = init
            ibins[i] = ibin
        end
    end
    binmass, ibins .- first(ibins)
end

function find_nearest_peak(::LocalMaxima, convolution, i, k)
    j = findfirst(>(0.5 * maximum(k)), k)
    ihwhm = floor(Int, length(k) / 2) - j + 1
    peak = [convolution[i], convolution[i]]
    dir = [true, true]
    start = [false, false]
    ibin = [i, i]
    for j in 0:ihwhm
        if first(dir)
            if convolution[i - j - 1] > first(peak)
                start[begin] = true
                peak[begin] = convolution[i - j - 1]
            else
                dir[begin] = false
                ibin[begin] = i - j
            end
        end
        if last(dir)
            if convolution[i + j + 1] > last(peak) 
                start[end] = true
                peak[end] = convolution[i + j + 1]
            else
                dir[end] = false
                ibin[end] = i + j
            end
        end
        dir'start > 0 || break 
    end
    r = @. (!)(dir) * start
    id = if !any(start)
        ibin[begin]
    elseif all(r) && peak[begin] < peak[end]
        ibin[end]
    elseif all(r)
        ibin[begin]
    elseif first(r)
        ibin[begin]
    elseif last(r)
        ibin[end]
    else
        nothing 
    end
    id
    # lm, _ = findmin(convolution[id - ihwhm : id])
    # rm, _ = findmin(convolution[id : id + ihwhm])
    # if lm > 0.75 * convolution[id] || rm > 0.75 * convolution[id]
    #     nothing 
    # else
    #     id
    # end
end

function find_nearest_peak(::FWHMMaxima, convolution, i, k)
    j = findfirst(>(0.5), k)
    ihwhm = floor(Int, length(k) / 2) - j + 1
    range = i - ihwhm : i + ihwhm
    _, v = findmax(convolution[range])
    mbin = range[v]
    n = ihwhm * 2
    if mbin == i - ihwhm
        while n > 0 && convolution[mbin] < convolution[mbin - 1]
            mbin -= 1
            n -= 1
            if n == 0
                mbin = nothing
                break 
            end
            if mbin == firstindex(convolution)
                break 
            end
        end
    elseif mbin == i + ihwhm
        while n > 0 && convolution[mbin] < convolution[mbin + 1]
            mbin += 1
            n -= 1
            if n == 0
                mbin = nothing
                break 
            end
            if mbin == lastindex(convolution)
                break 
            end
        end
    else
        lrange = mbin - ihwhm : mbin
        rrange = mbin : mbin + ihwhm
        lm, lmbin = findmin(convolution[lrange])
        rm, rmbin = findmin(convolution[rrange])
        lmbin = lrange[lmbin]
        rmbin = rrange[rmbin]
        if lm > 0.6 * convolution[mbin]
            while n > 0 && convolution[lmbin] < convolution[lmbin - 1]
                lmbin -= 1
                n -= 1
                if lmbin == firstindex(convolution)
                    break 
                end
            end
            if convolution[lmbin] < convolution[mbin]
                lmbin = -Inf
            end
        else
            lmbin = -Inf
        end
        if rm > 0.6 * convolution[mbin]
            while n > 0 && convolution[rmbin] < convolution[rmbin + 1]
                rmbin += 1
                n -= 1
                if rmbin == lastindex(convolution)
                    break 
                end
            end
            if convolution[rmbin] < convolution[mbin]
                rmbin = Inf
            end
        else
            rmbin = Inf
        end
        if isinf(rmbin) && isinf(lmbin)
            mbin = mbin
        elseif rmbin - mbin > mbin - lmbin
            mbin = lmbin 
        else
            mbin = rmbin
        end
    end
    mbin
end