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
    sp = string.(propertynames(mztable))
    colmz = findlastcol(sp, "MZ")
    colab = findlastcol(sp, "Abundance")
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
    colmz = findlastcol(string.(propertynames(mztable)), "MZ")
    if !isnothing(mz_range)
        lower_mz, upper_mz = mz_range 
        id = findall(x -> lower_mz < x < upper_mz, getproperty(mztable, colmz))
        mztable = mztable[id]
    end
    mztable
end
AllIons(mztable) = AllIons(nothing, mztable)
AllIons(mz_range, spec::Spectrum) = AllIons(mz_range, spec.table)

"""
    Isolation(msanalyzer, mztable; stage = nothing, threshold = rcrit(1e-4)) -> Table
    Isolation(msanalyzer, spectrum; stage = nothing, threshold = rcrit(1e-4)) -> Table

Isolating target ion(s) with specific m/z values and resolutions to enter the next MS stage. 

* `msanalyzer::AbstractMSAnalyzer`: a MS analyzer to perform MS filtering. See documentation of specific analyzer for detailed settings.
* `mztable::Table`: a table containing columns    
    * `ID`: ID tuples. Each elements represents ID number of ions of each MS stage. 
    * `Abundance1`, `Abundance2`, ..., `Abundancen`. 
    * `MZ1`, `MZ2`, ..., `MZn`. 
* `spectrum::Spectrum`: `spectrum.table` is utilized as `mztable`.
* `stage::Int`: MS stage. Default `nothing` for the last MS stage.
* `threshold` can be a number or criteria (absolute and/or relative to maximum), representing the lower limit of abundance. 
"""
function Isolation(msanalyzer::AbstractMSAnalyzer, mztable::Table; stage = nothing, threshold = rcrit(1e-4)) 
    isempty(mztable) && return mztable
    sp = string.(propertynames(mztable))
    if isnothing(stage)
        colmz = findlastcol(sp, "MZ")
        colab = findlastcol(sp, "Abundance")
    else
        colmz = findcol(sp, "MZ", stage)
        colab = findcol(sp, "Abundance", stage)
    end
    _Isolation(msanalyzer, mztable, colmz, colab, threshold)
end
Isolation(msanalyzer::AbstractMSAnalyzer, spec::Spectrum; stage = nothing, threshold = rcrit(1e-4)) = Isolation(msanalyzer, spec.table; stage, threshold)

function _Isolation(msanalyzer::AbstractMSAnalyzer, mztable::Table, colmz::Symbol, colab::Symbol, threshold)
    params = [window_parameter(msanalyzer, mz) for mz in vectorize(msanalyzer.mz)]
    ab = map(mztable) do r
        getproperty(r, colab) * maximum([msanalyzer.window(getproperty(r, colmz), param...) for param in params])
    end
    tab = minimum(makecrit_value(crit(threshold), maximum(ab)))
    id = findall(>(tab), ab)
    Table(mztable; [colab => ab]...)[id]
end

@deprecate TargetIon Isolation 

"""
    SelectedIonMonitor(transitiontable, mztable; threshold = rcrit(1e-4)) -> Table
    SelectedIonMonitor(transitiontable, spectrum; threshold = rcrit(1e-4)) -> Table

Selected ion monitoring. 

* `transitiontable::Table`: each row represents a transition. Use column `Transition` (optional) for specifying transition name. Other columns must be in analysis-fragmentation-analysis order. Analysis columns contain MS analyzers and Fragmentation columns contain producttables (See `Fragmentation` for detail).
* `mztable::Table`: a table containing columns    
    * `ID`: ID tuples. Each elements represents ID number of ions of each MS stage. 
    * `Abundance1`
    * `MZ1`
* `spectrum::Spectrum`: `spectrum.table` is utilized as `mztable`.
* `threshold` can be a number or criteria (absolute and/or relative to maximum), representing the lower limit of abundance. 
"""
function SelectedIonMonitor(transitiontable::Table, mztable::Table; threshold = rcrit(1e-4)) 
    isempty(mztable) && return transitiontable
    sp = string.(propertynames(mztable))
    colmz = findlastcol(sp, "MZ")
    colmz == :MZ1 || throw(ArgumentError("Require only single MZ column `MZ1."))
    colab = findlastcol(sp, "Abundance")
    colab == :Abundance1 || throw(ArgumentError("Require only single Abundance column `Abundance1."))
    transitions = Table(transitiontable; Transition = nothing)
    tables = _SelectedIonMonitor(columns(transitions), mztable, colmz, colab, threshold)
    Table(transitiontable; MZTable = tables)
end
SelectedIonMonitor(transitiontable::Table, spec::Spectrum; threshold = rcrit(1e-4)) = SelectedIonMonitor(transitiontable, spec.table; threshold)

function _SelectedIonMonitor(transitions, mztable::Table, colmz::Symbol, colab::Symbol, threshold)
    length(transitions) % 2 < 1 && throw(ArgumentError("Transitions must have 'analysis-fragmentation-analysis...` pattern."))
    isolation = true 
    mztable = repeat([mztable], length(first(transitions)))
    for transition in transitions
        if isolation
            mztable = [_Isolation(trans, mzt, colmz, colab, threshold) for (trans, mzt) in zip(transition, mztable)]
        else
            mztable = [Fragmentation(trans, mzt; threshold) for (trans, mzt) in zip(transition, mztable)]
            colmz = Symbol(string("MZ", parse(Int, string(match(r"\d+", string(colmz)).match)) + 1))
            colab = Symbol(string("Abundance", parse(Int, string(match(r"\d+", string(colab)).match)) + 1))
        end
        isolation = !isolation
    end
    mztable
end

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
    * `Abundance1`, `Abundance2`, ..., `Abundancen`. The last column will be utilized.
    * `MZ1`, `MZ2`, ..., `MZn`. The last column will be utilized.
* `spectrum::Spectrum`: `spectrum.table` is utilized as `precursortable`.
* `threshold` can be a number or criteria, representing the lower limit of abundance (absolute and/or relative to maximal value of each spectrum). 
"""
function Fragmentation(producttable::Table, mztable::Table; threshold = rcrit(1e-4))
    # group -> precursor_table, elements_precursor
    if !in(:ID, propertynames(producttable)) && in(:Chemical, propertynames(producttable))
        producttable = match_chemical(mztable, producttable; colexp = :Chemical, collib = :Chemical)
    elseif !in(:ID, propertynames(producttable))
        throw(ArgumentError("No column `ID` or `Chemical` in product_table."))
    end
    :ID in propertynames(mztable) || throw(ArgumentError("No column `ID` in precursor_table."))
    :Chemical in propertynames(mztable) || throw(ArgumentError("No column `Chemical` in precursor_table."))
    :Abundance1 in propertynames(mztable) || throw(ArgumentError("No column `Abundance1`, ..., `Abundacnen` in precursor_table."))
    :MZ1 in propertynames(mztable) || throw(ArgumentError("No column `MZ1`, ..., `MZn` in precursor_table."))
    :Product in propertynames(producttable) || throw(ArgumentError("No column `Product` in product_table."))
    allequal(msstage, mztable.Chemical) || throw(ArgumentError("Chemicals have to be in the same MS stage."))
    all(x -> all(y -> msstage(y) < 2, x), producttable.Product) || throw(ArgumentError("Products should not be MS/MS pairs."))
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
    peak_table(spectrum::Spectrum; alg = LocalMaxima(), abundance = 1, abtype = :max, threshold = rcrit(1e-4)) -> Table
    peak_table(transitiontable::Table; groupedisotopomers = true, isotope = "[13C]") -> Table

Extract peaks from a spectrum or SIM.

* `abundance` sets the abundance of the peak specified by `abtype`. 
* `abtype`
    * `:max`: the largest peak.
    * `:list`: sum of listed peaks.
    * `:raw`: no abundance normalization.
* `threshold` can be a number or criteria (absolute and/or relative to `abundance`), representing the lower limit of abundance. 
* `groupedisotopomers`: whether group isotopologues by isotopomer state based on `isotope`.
* `isotope::String`: minor isotope.
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
    c = length(first(mztable.Convolution)) ÷ 2 + 1
    gmztable = group(getproperty(:Max_bin), mztable)
    tuples = map(pairs(gmztable)) do (ibin, smztable) 
        cab = map(smztable) do r 
            delta = ibin - r.Bin_id
            r.Convolution[c + delta]
        end
        (; Chemical = Isobars(getproperty(smztable, :Chemical), cab), MZ = initial_mass + (ibin - 1) * binsize, Abundance = spectrum[ibin])
    end
    table = Table([t for t in tuples])
    normalize_abundance!(table.Abundance, abundance, abtype, [:max, :list, :raw])
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(table.Abundance)))
    table[findall(>=(abundance_cutoff), table.Abundance)]
end

function peak_table(transitiontable::Table; groupedisotopomers = true, isotope = "[13C]")
    :MZTable in propertynames(transitiontable) || throw(ArgumentError("No column`MZTable` in transitiontable"))
    mztables = filter(!isempty, transitiontable.MZTable)
    if groupedisotopomers
        tables = map(mztables) do table 
            group_isotopologues(table; isotope)
        end
    else
        tables = mztables
    end
    sp = string.(propertynames(first(tables)))
    colmz = findallcol(sp, "MZ")
    colab = findallcol(sp, "Abundance")
    gt = map(tables) do mztable 
        isobar = Isobars(mztable.Chemical, hcat(getproperty.(Ref(mztable), colab)...))
        transitions = chemicaltransition.(mztable.Chemical)
        if groupedisotopomers
            (; Chemical = isobar, [c => mean(getproperty(mztable, c), weights(getproperty(mztable, d))) for (c, d) in zip(colmz, colab)]..., [c => sum(getproperty(mztable, c)) for c in colab]...)
        else
            uid = [[findfirst(x -> x == t, c) for t in unique(c)] for c in zip(transitions...)]
            (; Chemical = isobar, [c => mean(getproperty(mztable, c)[i], weights(getproperty(mztable, d)[i])) for (i, c, d) in zip(uid, colmz, colab)]..., [c => sum(getproperty(mztable, c)[i]) for (i, c) in zip(uid, colab)]...)
        end
    end
    if :Transition in propertynames(transitiontable)
        Table(Table(; Transition = transitiontable.Transition), Table(collect(NamedTuple, gt)))
    else
        Table(collect(NamedTuple, gt))
    end
end

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

function find_nearest_peak(alg::LocalMaxima, convolution, i, k)
    j = findfirst(>(alg.threshold * maximum(k)), k)
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