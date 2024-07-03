function isotope_abundance(tbl::Table; colchemical = :Ion, colabundance = :Abundance, abundance = 1, colpreserve = setdiff(propertynames(tbl), [colchemical, colabundance]), kwargs...)
    colchemical in propertynames(tbl) || throw(ArgumentError("Column `$colchemical` does not exist."))
    mabundance = colabundance in propertynames(tbl) ? getproperty(tbl, colabundance) : vectorize(abundance, length(tbl))
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
        mapreduce(vcat, tbl, mabundance) do r, m
            x = isotope_abundance(getproperty(r, colchemical), m; kwargs...)
            Table(x; (map(colpreserve) do p
                p => vectorize(getproperty(r, p), length(x))
            end)...)
        end |> Table
    else
        t = Vector{Table}(undef, length(tbl))
        Threads.@threads for i in eachindex(t)
            x = isotope_abundance(getproperty(tbl[i], colchemical), mabundance[i]; kwargs...)
            t[i] = Table(x; (map(colpreserve) do p
                p => vectorize(getproperty(tbl[i], p), length(x))
            end)...)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end

function isotope_abundance(v::Vector, mabundance = 1; kwargs...)
    if length(v) < Threads.nthreads()
        mapreduce(vcat, v, vectorize(mabundance, length(v))) do x, m
            isotope_abundance(x, m; kwargs...)
        end |> Table
    else
        t = Vector{Table}(undef, length(v))
        mabundance = vectorize(mabundance, length(v))
        Threads.@threads for i in eachindex(t)
            t[i] = isotope_abundance(v[i], mabundance[i]; kwargs...)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end

getisotopes(ion::AbstractIon, mabundance = 1; abtype = :max, threshold = crit(mabundance * 1e-4, 1e-4), mz_tol = crit(0.01, 20e-6)) = 
    _isotope_abundance(ion, mabundance; abtype, threshold, mz_tol, table = false)
isotope_abundance(ion::AbstractIon, mabundance = 1; abtype = :max, threshold = crit(mabundance * 1e-4, 1e-4), mz_tol = crit(0.01, 20e-6)) = 
    _isotope_abundance(ion, mabundance; abtype, threshold, mz_tol, table = true)
function _isotope_abundance(ion::AbstractIon, mabundance = 1; abtype = :max, threshold = crit(mabundance * 1e-4, 1e-4), mz_tol = crit(0.01, 20e-6), table = true)
    it = _isotope_abundance(chemicalformula(ion), mabundance; abtype, threshold = threshold / 10)
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
    id = findall(>=(maximum(crit(threshold)(mabundance))), ab)
    table ? Table(; Ion = map((x, y) -> transform_ions(ion, x, y), e[id], a[id]), MZ = m[id], Abundance = ab[id]) : map((x, y) -> transform_ions(ion, x, y), e[id], a[id])
end

getisotopes(x::AbstractChemical, mabundance = 1; abtype = :max, threshold = crit(mabundance * 1e-4, 1e-4)) = 
    _isotope_abundance(x, mabundance; abtype, threshold, table = false)
isotope_abundance(x::AbstractChemical, mabundance = 1; abtype = :max, threshold = crit(mabundance * 1e-4, 1e-4)) = 
    _isotope_abundance(x, mabundance; abtype, threshold, table = true)
function _isotope_abundance(x::AbstractChemical, mabundance = 1; abtype = :max, threshold = crit(mabundance * 1e-4, 1e-4), table = true)
    it = _isotope_abundance(chemicalformula(x), mabundance; abtype, threshold)
    table ? Table(; Chemical = transform_chemicals(x, it.Element), MW = it.MW, Abundance = it.Abundance) : transform_chemicals(x, it.Element)
end

getisotopes(x::AbstractString, mabundance = 1; abtype = :max, threshold = crit(mabundance * 1e-4, 1e-4)) = map(chemicalformula, _isotope_abundance(x, mabundance; abtype, threshold, table = false))
function isotope_abundance(x::AbstractString, mabundance = 1; abtype = :max, threshold = crit(mabundance * 1e-4, 1e-4))
    it = _isotope_abundance(x, mabundance; abtype, threshold)
    Table(; Formula = map(chemicalformula, it.Element), MW = it.MW, Abundance = it.Abundance)    
end
function _isotope_abundance(x::AbstractString, mabundance = 1; abtype = :max, threshold = crit(mabundance * 1e-4, 1e-4), table = true)
    thresh = maximum(crit(threshold)(mabundance))
    el = parse_compound(x)
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
    ei = mapreduce(vcat, collect(keys(eld))) do i
        e = get(ELEMENTS, i, i)
        map(get(ISOTOPES, e, i)) do x
            (e, x)
        end
    end
    filter!(x -> !=(x...), ei)
    ab = abundance(el)
    threshn = @match abtype begin
        :all => thresh
        _    => thresh * ab
    end
    ite, ita = rec_isotopes!([vcat([k => v for (k, v) in eld], [k => v for (k, v) in itd])], [mabundance * ab], eld, itd, ei, 1, mabundance * ab, threshn)
    ita = abtype == :max ? ita ./ maximum(ita) .* mabundance : abtype == :mono ? ita ./ first(ita) .* mabundance : ita
    itm = map(mw, ite)
    id = sortperm(itm)
    filter!(x -> >=(ita[x], thresh), id)
    ite = ite[id]
    ita = ita[id]
    itm = itm[id]
    table ? Table(; Element = ite, MW = itm, Abundance = ita) : ite
end

abundance(x) = abundance(chemicalformula(x))
abundance(formula::AbstractString) = abundance(parse_compound(interpret_isotope(formula)))
function abundance(ite::Union{Vector, Dict})
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

function abundance(prev_abundance, old_element, new_element, nold, nnew, delta)
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
            ab = abundance(prev_abundance, e, i, el[e], new_it[i], 1)
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


