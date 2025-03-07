# Use Isotopes wrapper?
"""
    isotopologues(chemical::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mw_tol = crit(0.01, 20e-6))
    isotopologues(formula::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mw_tol = crit(0.01, 20e-6), net_charge = 0)

Isotopologues of a single `chemical` or `formula` (converted to `Chemical` by `parse_chemical`). 

* `abundance` sets the abundance of the isotope specified by `abtype`. 
    * `:max`: the most abundant isotopologue
    * `:input`: the input isotopologue
    * other: the final abundances are repressented as proportion.
* `threshold` can be a number or criteria, representing the lower limit of abundance. 
* `isobaric` determines whether groups isobars and creates `Isobars` or not.
* `mz_tol`: tolerance of m/z for isobars.
* `mw_tol`: tolerance of molecular weight for isobars.
* `net_charge`: charges (positive or negative) of `formula`.
"""
isotopologues(cc::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mw_tol = crit(0.01, 20e-6)) = 
    _isotopologues(cc, abundance; abtype, threshold, isobaric, mz_tol, mw_tol, table = false)

isotopologues(::Isobars, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mw_tol = crit(0.01, 20e-6)) = throw(ArgumentError("`Isobars` is not supported by `isotopologues`"))
isotopologues(::Isotopomers, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mw_tol = crit(0.01, 20e-6)) = throw(ArgumentError("`Isotopomers` is not supported by `isotopologues`"))
isotopologues(formula::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mw_tol = crit(0.01, 20e-6), net_charge = 0) = 
    _isotopologues(formula, abundance; abtype, threshold, isobaric, mz_tol, mw_tol, net_charge, table = false)

"""
    isotopologues_table(chemical::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mw_tol = crit(0.01, 20e-6))
    isotopologues_table(formula::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mw_tol = crit(0.01, 20e-6), net_charge = 0)
    isotopologues_table(tbl::Table; colchemical = :Chemical, colabundance = :Abundance, abundance = 1, colpreserve = setdiff(propertynames(tbl), [colchemical, colabundance]), kwargs...)
    isotopologues_table(v::Vector, abundance = 1; kwargs...)

A `Table` of isotopologues of single `chemical` or `formula` (converted to `Chemical` by `parse_chemical`), or multiple chemicals in `v` or column `colchemical` of `tbl`. 

* `abundance` sets the abundance of the isotope specified by `abtype`. 
    * `:max`: the most abundant isotopologue
    * `:input`: the input isotopologue
    * other: the final abundances are repressented as proportion.
* `threshold` can be a number or criteria, representing the lower limit of abundance. 
* `isobaric` determines whether groups isobars and creates `Isobars` or not.
* `mz_tol`: tolerance of m/z for isobars.
* `mw_tol`: tolerance of molecular weight for isobars.
* `net_charge`: charges (positive or negative) of `formula`.
"""
isotopologues_table(cc::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mw_tol = crit(0.01, 20e-6)) = 
    _isotopologues(cc, abundance; abtype, threshold, isobaric, mz_tol, mw_tol, table = true)

isotopologues_table(formula::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mw_tol = crit(0.01, 20e-6), net_charge = 0) = 
    _isotopologues(formula, abundance; abtype, threshold, isobaric, mz_tol, mw_tol, net_charge, table = true)

function _isotopologues(formula::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mw_tol = crit(0.01, 20e-6), net_charge = 0, table = true)
    if !isobaric
        it = _isotope_abundance(formula, abundance; abtype, threshold, net_charge)
        cs = map(chemicalformula, it.Element)
        ms = it.MW ./ max(1, abs(net_charge))
        table || return cs
        if charge(cc) == 0 
            tbl = Table(; Isotopologues = cs, MW = ms, Abundance = it.Abundance) 
        else
            tbl = Table(; Isotopologues = cs, MZ = ms, Abundance = it.Abundance) 
        end
        return tbl
    end
    it = _isotope_abundance(formula, abundance; abtype, threshold = threshold / 10, net_charge)
    # Table(; Formula = map(chemicalformula, it.Element), MW = it.MW, Abundance = it.Abundance)   
    m_tol = crit(real_interval(net_charge == 0 ? mw_tol : mz_tol))
    it.MW ./= max(1, abs(net_charge))
    e = typeof(it.Element)[]
    m = Float64[-Inf]
    a = Vector{Float64}[]
    for x in it
        if any(r -> in(x.MW - last(m), r), m_tol(last(m)))
            m[end] = (sum(last(a)) * m[end] + x.Abundance * x.MW) / (sum(last(a)) + x.Abundance)
            push!(last(e), x.Element)
            push!(last(a), x.Abundance)
        else
            push!(e, [x.Element])
            push!(a, [x.Abundance])
            push!(m, x.MW)
        end
    end
    popfirst!(m)
    for (x, y) in zip(e, a)
        ord = sortperm(y; rev = true)
        x .= x[ord]
    end
    ab = map(sum, a)
    id = findall(>=(maximum(crit(threshold)(abundance))), ab)
    if table 
        if net_charge == 0 
            Table(; Isotopologues = map(x -> chemicalformula.(x), e[id]), MW = m[id], Abundance = ab[id]) 
        else
            Table(; Isotopologues = map(x -> chemicalformula.(x), e[id]), MZ = m[id], Abundance = ab[id]) 
        end
    else 
        map(x -> chemicalformula.(x), e[id])
    end 
end

function _isotopologues(cc::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mw_tol = crit(0.01, 20e-6), table = true)
    # if charge(cc) == 0 
    #     it = _isotope_abundance(chemicalformula(cc), abundance; abtype, threshold)
    #     return table ? Table(; Chemical = Isotopomers.(x, it.Element), MW = it.MW, Abundance = it.Abundance) : Isotopomers.(x, it.Element)
    # end
    # thresh / 10, to make error of abundance less than 10%
    if !isobaric
        it = _isotope_abundance(chemicalformula(cc), abundance; abtype, threshold, net_charge = charge(cc))
        cs = Isotopomers.(cc, it.Element)
        ms = it.MW ./ max(1, ncharge(cc))
        table || return cs
        if charge(cc) == 0 
            tbl = Table(; Isotopologues = cs, MW = ms, Abundance = it.Abundance) 
        else
            tbl = Table(; Isotopologues = cs, MZ = ms, Abundance = it.Abundance) 
        end
        return tbl
    end
    it = _isotope_abundance(chemicalformula(cc), abundance; abtype, threshold = threshold / 10, net_charge = charge(cc))
    # it = Table(it; Element = map(x -> loss_elements!(unique_elements(x), adductelements(cc)), it.Element))
    m_tol = crit(real_interval(charge(cc) == 0 ? mw_tol : mz_tol))
    it.MW ./= max(1, ncharge(cc))
    e = typeof(it.Element)[]
    m = Float64[-Inf]
    a = Vector{Float64}[]
    for x in it
        if any(r -> in(x.MW - last(m), r), m_tol(last(m)))
            m[end] = (sum(last(a)) * m[end] + x.Abundance * x.MW) / (sum(last(a)) + x.Abundance)
            push!(last(e), x.Element)
            push!(last(a), x.Abundance)
        else
            push!(e, [x.Element])
            push!(a, [x.Abundance])
            push!(m, x.MW)
        end
    end
    popfirst!(m)
    for (x, y) in zip(e, a)
        ord = sortperm(y; rev = true)
        x .= x[ord]
    end
    ab = map(sum, a)
    id = findall(>=(maximum(crit(threshold)(abundance))), ab)
    if table 
        if charge(cc) == 0 
            Table(; Isotopologues = map((x, y) -> Isobars(Isotopomers.(cc, x), y), e[id], a[id]), MW = m[id], Abundance = ab[id]) 
        else
            Table(; Isotopologues = map((x, y) -> Isobars(Isotopomers.(cc, x), y), e[id], a[id]), MZ = m[id], Abundance = ab[id]) 
        end
    else 
        map((x, y) -> Isobars(Isotopomers.(cc, x), y), e[id], a[id])
    end
end

# loewest level
function _isotope_abundance(x::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), net_charge = 0, table = true)
    thresh = maximum(crit(threshold)(abundance))
    el = chemicalelements(x)
    # record elements change
    eld = Dictionary{String, Int}()
    itd = Dictionary{String, Int}()
    for (e, n) in el
        if get(ELEMENTS[:PARENTS], e, e) == e
            get!(eld, e, 0)
            eld[e] += n
        else
            get!(itd, e, 0)
            itd[e] += n
        end
    end
    # element => isoptope pairs
    ei = mapreduce(vcat, collect(keys(eld))) do i
        e = get(ELEMENTS[:PARENTS], i, i)
        map(get(ELEMENTS[:ISOTOPES], e, i)) do x
            (e, x)
        end
    end
    filter!(x -> !=(x...), ei)
    sort!(ei; by = x -> ELEMENTS[:ABUNDANCE][last(x)], rev = true)
    ab = isotopicabundance(el; ignore_isotopes = true)
    threshn = @match abtype begin
        :input => thresh * ab
        :max  => thresh * ab
        _     => thresh
    end
    # serve abundance as sums
    ite, ita = rec_isotopes!([unique_elements(el)], [abundance * ab], eld, itd, ei, 1, abundance * ab, threshn)
    ita = abtype == :max ? ita ./ maximum(ita) .* abundance : abtype == :input ? ita ./ first(ita) .* abundance : ita
    itm = map(mw, ite, repeat([net_charge], length(ite)))
    id = sortperm(itm)
    filter!(x -> >=(ita[x], thresh), id)
    ite = ite[id]
    ita = ita[id]
    itm = itm[id]
    table ? Table(; Element = ite, MW = itm, Abundance = ita) : ite
end

function isotopologues_table(tbl::Table; colchemical = :Chemical, colabundance = :Abundance, abundance = 1, colpreserve = setdiff(propertynames(tbl), [colchemical, colabundance]), kwargs...)
    colchemical in propertynames(tbl) || throw(ArgumentError("Column `$colchemical` does not exist."))
    abundance = colabundance in propertynames(tbl) ? getproperty(tbl, colabundance) : vectorize(abundance, length(tbl))
    del = Int[]
    for (i, p) in enumerate(colpreserve)
        if !in(p, propertynames(tbl))
            @warn "Column `$p` does not exist. Ignore this column."
            push!(del, i)
        elseif p == :Chemical || p == :MZ || p == :Abundance
            @warn "Column `$p` is preserved. Ignore this column."
            push!(del, i)
        end
    end
    colpreserve = deleteat!(collect(colpreserve), del)
    if length(tbl) < Threads.nthreads()
        mapreduce(vcat, tbl, abundance) do r, m
            x = isotopologues_table(getproperty(r, colchemical), m; kwargs...)
            Table(x; (map(colpreserve) do p
                p => vectorize(getproperty(r, p), length(x))
            end)...)
        end |> Table
    else
        t = Vector{Table}(undef, length(tbl))
        Threads.@threads for i in eachindex(t)
            x = isotopologues_table(getproperty(tbl[i], colchemical), abundance[i]; kwargs...)
            t[i] = Table(x; (map(colpreserve) do p
                p => vectorize(getproperty(tbl[i], p), length(x))
            end)...)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end

function isotopologues_table(v::Vector, abundance = 1; kwargs...)
    if length(v) < Threads.nthreads()
        mapreduce(vcat, v, vectorize(abundance, length(v))) do x, m
            isotopologues_table(x, m; kwargs...)
        end |> Table
    else
        t = Vector{Table}(undef, length(v))
        abundance = vectorize(abundance, length(v))
        Threads.@threads for i in eachindex(t)
            t[i] = isotopologues_table(v[i], abundance[i]; kwargs...)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end

isotopologues_table(::Isobars, abundance = 1; kwargs...) = throw(ArgumentError("`Isobars` is not supported by `isotopologues_table`"))
isotopologues_table(::Isotopomers, abundance = 1; kwargs...) = throw(ArgumentError("`Isotopomers` is not supported by `isotopologues_table`"))

"""
    isotopicabundance(chemical::AbstractChemical; ignore_isotopes = false)
    isotopicabundance(formula::AbstractString; ignore_isotopes = false)
    isotopicabundance(elements::Union{<: Vector, Dictionary}; ignore_isotopes = false)

Compute isotopic abundance of `chemical`, `formula`, vector of element-number pairs or dictionary mapping element to number. 
"""
isotopicabundance(cc::AbstractChemical; ignore_isotopes = false) = isotopicabundance(chemicalformula(cc); ignore_isotopes)
isotopicabundance(formula::AbstractString; ignore_isotopes = false) = isotopicabundance(chemicalelements(formula); ignore_isotopes)
isotopicabundance(elements::Dictionary; ignore_isotopes = false) = isotopicabundance(pairs(elements); ignore_isotopes)
function isotopicabundance(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}; ignore_isotopes = false)
    el = Dictionary{String, Vector{Int}}()
    for (e, n) in elements
        haskey(ELEMENTS[:PARENTS], e) || continue
        k = ELEMENTS[:PARENTS][e]
        haskey(ELEMENTS[:ISOTOPES], k) || continue
        get!(el, k, zeros(Int, length(ELEMENTS[:ISOTOPES][k])))
        (!ignore_isotopes || k == e) && (el[k][findfirst(==(e), ELEMENTS[:ISOTOPES][k])] += n)
    end # use two vector?
    ab = 1
    for (e, v) in pairs(el)
        as = get.(Ref(ELEMENTS[:ABUNDANCE]), ELEMENTS[:ISOTOPES][e], 1)
        any(==(1), as) && continue
        n = sum(v)
        id = sortperm(v; rev = true)
        v = v[id]
        as = as[id]
        n1 = popfirst!(v)
        f = try 
            factorial(n, n1)
        catch 
            factorial(big(n), n1)
        end
        ab *= f * popfirst!(as) ^ n1
        isempty(as) && continue
        ab *= mapreduce(*, v, as) do x, y
            y ^ x / factorial(x)
        end
    end
    Float64(ab)
end

function update_abundance(prev_abundance, old_element, new_element, nold, nnew, delta)
    x = get(ELEMENTS[:ABUNDANCE], old_element, 1)
    y = get(ELEMENTS[:ABUNDANCE], new_element, 1)
    (x == 1 || y == 1) && return prev_abundance
    prev_abundance * factorial(nold, nold - delta) / factorial(nnew + delta, nnew) * (y / x) ^ delta
end

function rec_isotopes!(ite::Vector, ita::Vector, el::Dictionary, it::Dictionary, ei::Vector, ip::Int, prev_abundance, threshold)
    new_ip = -1
    for (e, i) in @views ei[ip:end]
        ne = get(el, e, 0)
        new_ip += 1
        if ne > 0
            new_el = deepcopy(el)
            new_it = deepcopy(it)
            new_el[e] -= 1
            get!(new_it, i, 0)
            ab = update_abundance(prev_abundance, e, i, el[e], new_it[i], 1)
            new_it[i] += 1
            ab < threshold && break
            push!(ite, add_elements(new_el, new_it))
            push!(ita, ab)
            rec_isotopes!(ite, ita, new_el, new_it, ei, ip + new_ip, ab, threshold)
        end
    end
    ite, ita
end