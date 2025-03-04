"""
    getisotopes(ion::AbstractIon, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), mz_tol = crit(0.01, 20e-6))
    getisotopes(m::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4))
    getisotopes(x::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4))

Isotopes of `ion`, `m` or `x` (converted to `Chemical` by `parse_chemical`). 

`abundance` sets the abundance of the isotope specified by `abtype`. 
* `:max`: the most abundant isotope
* `:mono`: the first isotope
* other: the final abundances are repressented as proportion.

`threshold` can be a number or criteria, representing the lower limit of abundance. 
"""
getisotopes(ion::AbstractIon, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), mz_tol = crit(0.01, 20e-6)) = 
    _isotope_abundance(ion, abundance; abtype, threshold, mz_tol, table = false)

getisotopes(m::Isobars, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), mz_tol = crit(0.01, 20e-6)) = throw(ArgumentError("`Isobars` is not supported by `getisotopes`"))

"""
    isotopetable(ion::AbstractIon, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), mz_tol = crit(0.01, 20e-6))
    isotopetable(m::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4))
    isotopetable(x::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4))
    isotopetable(tbl::Table; colchemical = :Ion, colabundance = :Abundance, abundance = 1, colpreserve = setdiff(propertynames(tbl), [colchemical, colabundance]), kwargs...)
    isotopetable(v::Vector, abundance = 1; kwargs...)
    isotopetable(v::Isobars, abundance = 1; kwargs...)

A `Table` of isotopes of single `ion`, `m` or `x` (converted to `Chemical` by `parse_chemical`), or multiple chemicals in `v` or column `colchemical` of `tbl`. 

`abundance` sets the abundance of the specific isotope determined by `abtype`.
* `:max`: the most abundant isotope
* `:mono`: the first isotope
* other: the final abundances are repressented as proportion.

Keyword argument `colpreserve` specifies columns to be presserved in new table. 

`threshold` can be a number or criteria, representing the lower limit of abundance. 
"""
isotopetable(ion::AbstractIon, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), mz_tol = crit(0.01, 20e-6)) = 
    _isotope_abundance(ion, abundance; abtype, threshold, mz_tol, table = true)
function _isotope_abundance(ion::AbstractIon, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), mz_tol = crit(0.01, 20e-6), table = true)
    # thresh / 10, to make error of abundance less than 10%
    it = _isotope_abundance(chemicalformula(ion), abundance; abtype, threshold = threshold / 10)
    it = Table(it; Element = map(x -> loss_elements!(unique_elements(x), adductelement(ion)), it.Element))
    mz_tol = crit(real_interval(mz_tol))
    it.MW ./= charge(ion)
    e = typeof(it.Element)[]
    m = Float64[-Inf]
    a = Vector{Float64}[]
    for x in it
        if any(r -> in(x.MW - last(m), r), mz_tol(last(m)))
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
    table ? Table(; Ion = map((x, y) -> Isobars(ionvariants(ion, x), y), e[id], a[id]), MZ = m[id], Abundance = ab[id]) : map((x, y) -> Isobars(ionvariants(ion, x), y), e[id], a[id])
end

getisotopes(x::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4)) = 
    _isotope_abundance(x, abundance; abtype, threshold, table = false)
isotopetable(x::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4)) = 
    _isotope_abundance(x, abundance; abtype, threshold, table = true)
function _isotope_abundance(x::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), table = true)
    it = _isotope_abundance(chemicalformula(x), abundance; abtype, threshold)
    table ? Table(; Chemical = chemicalvariants(x, it.Element), MW = it.MW, Abundance = it.Abundance) : chemicalvariants(x, it.Element)
end

getisotopes(x::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4)) = map(chemicalformula, _isotope_abundance(x, abundance; abtype, threshold, table = false))
function isotopetable(x::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4))
    it = _isotope_abundance(x, abundance; abtype, threshold)
    Table(; Formula = map(chemicalformula, it.Element), MW = it.MW, Abundance = it.Abundance)    
end
function _isotope_abundance(x::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), table = true)
    thresh = maximum(crit(threshold)(abundance))
    el = parse_compound(x)
    # record elements change
    eld = Dict{String, Int}()
    itd = Dict{String, Int}()
    for (e, n) in el
        if get(ELEMENTS, e, e) == e
            get!(eld, e, 0)
            eld[e] += n
        else
            get!(itd, e, 0)
            itd[e] += n
        end
    end
    # element => isoptope pairs
    ei = mapreduce(vcat, collect(keys(eld))) do i
        e = get(ELEMENTS, i, i)
        map(get(ISOTOPES, e, i)) do x
            (e, x)
        end
    end
    filter!(x -> !=(x...), ei)
    ab = isotopicabundance(el)
    threshn = @match abtype begin
        :mono => thresh * ab
        :max  => thresh * ab
        _     => thresh
    end
    # serve abundance as sums
    ite, ita = rec_isotopes!([vcat([k => v for (k, v) in eld], [k => v for (k, v) in itd])], [abundance * ab], eld, itd, ei, 1, abundance * ab, threshn)
    ita = abtype == :max ? ita ./ maximum(ita) .* abundance : abtype == :mono ? ita ./ first(ita) .* abundance : ita
    itm = map(mw, ite)
    id = sortperm(itm)
    filter!(x -> >=(ita[x], thresh), id)
    ite = ite[id]
    ita = ita[id]
    itm = itm[id]
    table ? Table(; Element = ite, MW = itm, Abundance = ita) : ite
end

function isotopetable(tbl::Table; colchemical = :Ion, colabundance = :Abundance, abundance = 1, colpreserve = setdiff(propertynames(tbl), [colchemical, colabundance]), kwargs...)
    colchemical in propertynames(tbl) || throw(ArgumentError("Column `$colchemical` does not exist."))
    abundance = colabundance in propertynames(tbl) ? getproperty(tbl, colabundance) : vectorize(abundance, length(tbl))
    del = Int[]
    for (i, p) in enumerate(colpreserve)
        if !in(p, propertynames(tbl))
            @warn "Column `$p` does not exist. Ignore this column."
            push!(del, i)
        elseif p == :Ion || p == :MZ || p == :Abundance
            @warn "Column `$p` is preserved. Ignore this column."
            push!(del, i)
        end
    end
    colpreserve = deleteat!(collect(colpreserve), del)
    if length(tbl) < Threads.nthreads()
        mapreduce(vcat, tbl, abundance) do r, m
            x = isotopetable(getproperty(r, colchemical), m; kwargs...)
            Table(x; (map(colpreserve) do p
                p => vectorize(getproperty(r, p), length(x))
            end)...)
        end |> Table
    else
        t = Vector{Table}(undef, length(tbl))
        Threads.@threads for i in eachindex(t)
            x = isotopetable(getproperty(tbl[i], colchemical), abundance[i]; kwargs...)
            t[i] = Table(x; (map(colpreserve) do p
                p => vectorize(getproperty(tbl[i], p), length(x))
            end)...)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end

function isotopetable(v::Vector, abundance = 1; kwargs...)
    if length(v) < Threads.nthreads()
        mapreduce(vcat, v, vectorize(abundance, length(v))) do x, m
            isotopetable(x, m; kwargs...)
        end |> Table
    else
        t = Vector{Table}(undef, length(v))
        abundance = vectorize(abundance, length(v))
        Threads.@threads for i in eachindex(t)
            t[i] = isotopetable(v[i], abundance[i]; kwargs...)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end

isotopetable(m::Isobars, abundance = 1; kwargs...) = isotopetable(n.chemicals, abundance; kwargs...)


"""
    isotopicabundance(m)
    isotopicabundance(formula::AbstractString)
    isotopicabundance(ite::Union{Vector, Dict})

Compute isotopic abundance of chemical `m`, formula, vector of element-number pairs or dictionary mapping element to number. 
"""
isotopicabundance(x) = isotopicabundance(chemicalformula(x))
isotopicabundance(formula::AbstractString) = isotopicabundance(parse_compound(transform_isotope_repr(formula)))
function isotopicabundance(ite::Union{Vector, Dict})
    el = Dict{String, Vector{Int}}()
    for (e, n) in ite
        haskey(ELEMENTS, e) || continue
        k = ELEMENTS[e]
        haskey(ISOTOPES, k) || continue
        get!(el, k, zeros(Int, length(ISOTOPES[k])))
        el[k][findfirst(==(e), ISOTOPES[k])] += n
    end # use two vector?
    ab = 1
    for (e, v) in el
        as = get.(Ref(ABUNDANCE), ISOTOPES[e], 1)
        any(==(1), as) && continue
        n = sum(v)
        n1 = popfirst!(v)
        ab *= factorial(n, n1) * popfirst!(as) ^ n1
        isempty(as) && continue
        ab *= mapreduce(*, v, as) do x, y
            y ^ x / factorial(x)
        end
    end
    ab
end

function update_abundance(prev_abundance, old_element, new_element, nold, nnew, delta)
    x = get(ABUNDANCE, old_element, 1)
    y = get(ABUNDANCE, new_element, 1)
    (x == 1 || y == 1) && return prev_abundance
    prev_abundance * factorial(nold, nold - delta) / factorial(nnew + delta, nnew) * (y / x) ^ delta
end

function rec_isotopes!(ite::Vector, ita::Vector, el::Dict, it::Dict, ei::Vector, ip::Int, prev_abundance, threshold)
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
            new = vcat([k => v for (k, v) in new_el], [k => v for (k, v) in new_it])
            push!(ite, new)
            push!(ita, ab)
            rec_isotopes!(ite, ita, new_el, new_it, ei, ip + new_ip, ab, threshold)
        end
    end
    ite, ita
end


