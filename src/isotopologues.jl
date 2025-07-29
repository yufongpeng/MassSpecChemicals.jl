# Use Isotopes wrapper?
"""
    isotopologues(chemical::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6))
    isotopologues(formula::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6), net_charge = 0)
    isotopologues(chemicalpair::ChemicalPair, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6))) 
    isotopologues(formula::Pair{<: AbstractString, AbstractString}, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), net_charge = (1, 1))

Isotopologues of a single `chemical` or `formula` (converted to `Chemical` by `parse_chemical`). 

* `abundance` sets the abundance of the isotope specified by `abtype`. 
    * `:max`: the most abundant isotopologue
    * `:input`: the input isotopologue
    * other: the final abundances are repressented as proportion.
* `threshold` can be a number or criteria, representing the lower limit of abundance. 
* `isobaric` determines whether groups isobars and creates `Isobars` or not.
* `mz_tol`: tolerance of m/z for isobars.
* `mm_tol`: tolerance of molecular mass for isobars.
* `net_charge`: charges (positive or negative) of `formula`.
"""
isotopologues(cc::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6)) = 
    _isotopologues(cc, abundance; abtype, threshold, isobaric, mz_tol, mm_tol, table = false)

isotopologues(::Isobars, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6)) = throw(ArgumentError("`Isobars` is not supported by `isotopologues`"))
isotopologues(::Isotopomers, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6)) = throw(ArgumentError("`Isotopomers` is not supported by `isotopologues`"))
isotopologues(formula::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6), net_charge = 0) = 
    _isotopologues(formula, abundance; abtype, threshold, isobaric, mz_tol, mm_tol, net_charge, table = false)

isotopologues(cc::ChemicalPair, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6))) = 
    _isotopologues(cc, abundance; abtype, threshold, isobaric, mz_tol, mm_tol, table = false)

isotopologues(formula::Pair{<: AbstractString, AbstractString}, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), net_charge = (1, 1)) = 
    _isotopologues(formula, abundance; abtype, threshold, isobaric, mz_tol, mm_tol, net_charge, table = false)
# isotopolgues_msn
"""
    isotopologues_table(chemical::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6))
    isotopologues_table(formula::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6), net_charge = 0)
    isotopologues_table(chemicalpair::ChemicalPair, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6))) 
    isotopologues_table(formula::Pair{<: AbstractString, <: AbstractString}, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), net_charge = (1, 1)) 
    isotopologues_table(tbl::Table; colchemical = :Chemical, colabundance = :Abundance, abundance = 1, colpreserve = setdiff(propertynames(tbl), [colchemical, colabundance]), kwargs...)
    isotopologues_table(v::Vector, abundance = 1; kwargs...)

A `Table` of isotopologues of single `chemical` or `formula` (converted to `Chemical` by `parse_chemical`), chemical pairs, or multiple chemicals in `v` or column `colchemical` of `tbl`. 

* `abundance` sets the abundance of the isotope specified by `abtype`. 
    * `:max`: the most abundant isotopologue
    * `:input`: the input isotopologue
    * other: the final abundances are repressented as proportion.
* `threshold` can be a number or criteria, representing the lower limit of abundance. 
* `isobaric` determines whether groups isobars and creates `Isobars` or not.
* `mz_tol`: tolerance of m/z for isobars.
* `mm_tol`: tolerance of molecular mass for isobars.
* `net_charge`: charges (positive or negative) of `formula`.
"""
isotopologues_table(cc::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6)) = 
    _isotopologues(cc, abundance; abtype, threshold, isobaric, mz_tol, mm_tol, table = true)

isotopologues_table(formula::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6), net_charge = 0) = 
    _isotopologues(formula, abundance; abtype, threshold, isobaric, mz_tol, mm_tol, net_charge, table = true)
    
isotopologues_table(cc::ChemicalPair, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6))) = 
    _isotopologues(cc, abundance; abtype, threshold, isobaric, mz_tol, mm_tol, table = true)

isotopologues_table(formula::Pair{<: AbstractString, <: AbstractString}, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), net_charge = (1, 1)) = 
    _isotopologues(formula, abundance; abtype, threshold, isobaric, mz_tol, mm_tol, net_charge, table = true)

function _isotopologues(formula::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6), net_charge = 0, colision = 0, table = true)
    if !isobaric
        it = _isotope_abundance(formula, abundance; abtype, threshold, net_charge)
        ms = it.Mass ./ max(1, abs(net_charge))
        colision > 0 && return (it.Element, ms, it.Abundance)
        cs = map(chemicalformula, it.Element)
        table || return cs
        if net_charge == 0 
            tbl = Table(; Isotopologues = cs, Mass = ms, Abundance = it.Abundance) 
        else
            tbl = Table(; Isotopologues = cs, MZ = ms, Abundance = it.Abundance) 
        end
        return tbl
    end
    it = _isotope_abundance(formula, abundance; abtype, threshold = threshold / 10, net_charge)
    # Table(; Formula = map(chemicalformula, it.Element), Mass = it.Mass, Abundance = it.Abundance)   
    m_tol = crit(net_charge == 0 ? mm_tol : mz_tol)
    it.Mass ./= max(1, abs(net_charge))
    e, m, a = isobaric_sum(it.Element, it.Mass, it.Abundance, m_tol)
    ab = map(sum, a)
    id = findall(>=(maximum(makecrit_value(crit(threshold), abundance))), ab)
    if colision > 0
        (e[id], m[id], a[id])
    elseif table 
        if net_charge == 0 
            Table(; Isotopologues = map(x -> chemicalformula.(x), e[id]), Mass = m[id], Abundance = ab[id]) 
        else
            Table(; Isotopologues = map(x -> chemicalformula.(x), e[id]), MZ = m[id], Abundance = ab[id]) 
        end
    else 
        map(x -> chemicalformula.(x), e[id])
    end 
end

function _isotopologues(cc::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6), colision = 0, table = true)
    # if charge(cc) == 0 
    #     it = _isotope_abundance(chemicalformula(cc), abundance; abtype, threshold)
    #     return table ? Table(; Chemical = Isotopomers.(x, it.Element), Mass = it.Mass, Abundance = it.Abundance) : Isotopomers.(x, it.Element)
    # end
    # thresh / 10, to make error of abundance less than 10%
    if !isobaric
        it = _isotope_abundance(chemicalformula(cc), abundance; abtype, threshold, net_charge = charge(cc))
        ms = it.Mass ./ max(1, ncharge(cc))
        colision > 0 && return (it.Element, ms, it.Abundance)
        cs = Isotopomers.(cc, it.Element)
        table || return cs
        if charge(cc) == 0 
            tbl = Table(; Isotopologues = cs, Mass = ms, Abundance = it.Abundance) 
        else
            tbl = Table(; Isotopologues = cs, MZ = ms, Abundance = it.Abundance) 
        end
        return tbl
    end
    it = _isotope_abundance(chemicalformula(cc), abundance; abtype, threshold = threshold / 10, net_charge = charge(cc))
    # it = Table(it; Element = map(x -> loss_elements!(unique_elements(x), adductelements(cc)), it.Element))
    m_tol = crit(charge(cc) == 0 ? mm_tol : mz_tol)
    it.Mass ./= max(1, ncharge(cc))
    e, m, a = isobaric_sum(it.Element, it.Mass, it.Abundance, m_tol)
    ab = map(sum, a)
    id = findall(>=(maximum(makecrit_value(crit(threshold), abundance))), ab)
    if colision > 0
        (e[id], m[id], a[id])
    elseif table 
        if charge(cc) == 0 
            Table(; Isotopologues = map((x, y) -> Isobars(Isotopomers.(cc, x), y), e[id], a[id]), Mass = m[id], Abundance = ab[id]) 
        else
            Table(; Isotopologues = map((x, y) -> Isobars(Isotopomers.(cc, x), y), e[id], a[id]), MZ = m[id], Abundance = ab[id]) 
        end
    else 
        map((x, y) -> Isobars(Isotopomers.(cc, x), y), e[id], a[id])
    end
end

function _isotopologues(formula::Pair{<: AbstractString, <: AbstractString}, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), net_charge = (1, 1), colision = 1, table = true)
    any(==(0), net_charge) && any(!=(0), net_charge) && throw(ArgumentError("Charges must be all non-zero or all zero."))
    e, m, a = _isotopologues(first(formula), abundance; abtype, threshold, isobaric, mz_tol = first(mz_tol), mm_tol = first(mm_tol), net_charge = first(net_charge), colision, table)
    ef = chemicalelements(last(formula))
    if !isobaric
        ee = Pair[]
        mm = Float64[]
        mm1 = Float64[]
        aa = Float64[]
        for (ei, mi, ai) in zip(e, m, a)
            ek, mk, ak = isotopes_proportion(ei, unique_elements(ef), loss_elements(chemicalelements(first(formula)), ef), last(net_charge))
            id = sortperm(mk)
            ee = vcat(ee, Ref(ei) .=> ek[id])
            mm = vcat(mm, mk[id])
            mm1 = vcat(mm1, repeat([mi], length(ek)))
            aa = vcat(aa, ai .* ak[id])
        end
        colision -= 1
        colision > 0 && return (ee, mm1 .=> mm, aa)
        cs = map(x -> chemicalformula(first(x)) => chemicalformula(last(x)), ee)
        table || return cs
        if last(net_charge) == 0 
            tbl = Table(; Isotopologues = cs, Mass1 = mm1, Mass2 = mm, Abundance = aa) 
        else
            tbl = Table(; Isotopologues = cs, MZ1 = mm1, MZ2 = mm, Abundance = aa) 
        end
        return tbl
    end
    m_tol = crit(last(net_charge) == 0 ? last(mm_tol) : last(mz_tol))
    eee = Pair[]
    mm1 = Float64[]
    mmm = Float64[]
    aaa = Float64[]
    for (ei, mi, ai) in zip(e, m, a)
        ee = Pair[]
        mm = Float64[]
        aa = Float64[]
        for (ej, aj) in zip(ei, ai)
            ek, mk, ak = isotopes_proportion(ej, unique_elements(ef), loss_elements(chemicalelements(first(formula)), ef), last(net_charge))
            ee = vcat(ee, Ref(ej) .=> ek)
            mm = vcat(mm, mk)
            aa = vcat(aa, aj .* ak)
        end
        id = sortperm(mm)
        ee = ee[id]
        mm = mm[id]
        aa = aa[id]
        mm ./= max(1, abs(last(net_charge)))
        ee, mm, aa = isobaric_sum(ee, mm, aa, m_tol)
        eee = vcat(eee, ee)
        mmm = vcat(mmm, mm)
        mm1 = vcat(mm1, repeat([mi], length(ee)))
        aaa = vcat(aaa, map(sum, aa))
    end
    id = findall(>=(maximum(makecrit_value(crit(threshold), abundance))), aaa)
    colision -= 1
    if colision > 0
        (eee[id], mmm[id], aaa[id])
    elseif table 
        if last(net_charge) == 0 
            Table(; Isotopologues = map(x -> map(y -> Pair(chemicalformula(first(y)), chemicalformula(last(y))), x), eee[id]), Mass1 = mm1[id], Mass2 = mmm[id], Abundance = aaa[id]) 
        else
            Table(; Isotopologues = map(x -> map(y -> Pair(chemicalformula(first(y)), chemicalformula(last(y))), x), eee[id]), MZ1 = mm1[id], MZ2 = mmm[id], Abundance = aaa[id]) 
        end
    else 
        map(x -> map(y -> Pair(chemicalformula(first(y)), chemicalformula(last(y))), x), eee[id])
    end 
end

function _isotopologues(cc::ChemicalPair, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), colision = 1, table = true)
    cs = (ncharge(cc.precursor), ncharge(cc.product))
    any(==(0), cs) && any(!=(0), cs) && throw(ArgumentError("Charges must be all non-zero or all zero."))
    e, m, a = _isotopologues(cc.precursor, abundance; abtype, threshold, isobaric, mz_tol = first(mz_tol), mm_tol = first(mm_tol), colision, table)
    ef = chemicalelements(cc.product)
    if !isobaric
        ee = Pair[]
        mm = Float64[]
        mm1 = Float64[]
        aa = Float64[]
        for (ei, mi, ai) in zip(e, m, a)
            ek, mk, ak = isotopes_proportion(ei, unique_elements(ef), loss_elements(chemicalelements(cc.precursor), ef), last(cs))
            id = sortperm(mk)
            ee = vcat(ee, Ref(ei) .=> ek[id])
            mm = vcat(mm, mk[id])
            mm1 = vcat(mm1, repeat([mi], length(ek)))
            aa = vcat(aa, ai .* ak[id])
        end
        colision -= 1
        colision > 0 && return (ee, mm1 .=> mm, aa)
        cs = ChemicalPair.(Isotopomers.(cc.precursor, first.(ee)), Isotopomers.(cc.product, last.(ee)))
        table || return cs
        if last(cs) == 0 
            tbl = Table(; Isotopologues = cs, Mass1 = mm1, Mass2 = mm, Abundance = aa) 
        else
            tbl = Table(; Isotopologues = cs, MZ1 = mm1, MZ2 = mm, Abundance = aa) 
        end
        return tbl
    end
    m_tol = crit(last(cs) == 0 ? last(mm_tol) : last(mz_tol))
    eee = Pair[]
    mm1 = Float64[]
    mmm = Float64[]
    aaa = Vector{Float64}[]
    for (ei, mi, ai) in zip(e, m, a)
        ee = Pair[]
        mm = Float64[]
        aa = Float64[]
        for (ej, aj) in zip(ei, ai)
            ek, mk, ak = isotopes_proportion(ej, unique_elements(ef), loss_elements(chemicalelements(cc.precursor), ef), last(cs))
            ee = vcat(ee, Ref(ej) .=> ek)
            mm = vcat(mm, mk)
            aa = vcat(aa, aj .* ak)
        end
        id = sortperm(mm)
        ee = ee[id]
        mm = mm[id]
        aa = aa[id]
        mm ./= max(1, abs(last(cs)))
        ee, mm, aa = isobaric_sum(ee, mm, aa, m_tol)
        eee = vcat(eee, ee)
        mmm = vcat(mmm, mm)
        mm1 = vcat(mm1, repeat([mi], length(ee)))
        aaa = vcat(aaa, aa)
    end
    aba = map(sum, aaa)
    id = findall(>=(maximum(makecrit_value(crit(threshold), abundance))), aba)
    colision -= 1
    # Isotopologues = map((x, y) -> Isobars(Isotopomers.(cc, x), y), e[id], a[id])
    if colision > 0
        (eee[id], mmm[id], aba[id])
    elseif table 
        if last(cs) == 0 
            Table(; Isotopologues = map((x, a) -> Isobars(map(y -> ChemicalPair(Isotopomers(cc.precursor, first(y)), Isotopomers(cc.product, last(y))), x), a), eee[id], aaa[id]), Mass1 = mm1[id], Mass2 = mmm[id], Abundance = aba[id]) 
        else
            Table(; Isotopologues = map((x, a) -> Isobars(map(y -> ChemicalPair(Isotopomers(cc.precursor, first(y)), Isotopomers(cc.product, last(y))), x), a), eee[id], aaa[id]), MZ1 = mm1[id], MZ2 = mmm[id], Abundance = aba[id]) 
        end
    else 
        map((x, a) -> Isobars(map(y -> ChemicalPair(Isotopomers(cc.precursor, first(y)), Isotopomers(cc.product, last(y))), x), a), eee[id], aaa[id])
    end 
end

# loewest level
function _isotope_abundance(x::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), net_charge = 0, table = true)
    thresh = maximum(makecrit_value(crit(threshold), abundance))
    el = chemicalelements(x)
    # record elements change
    eld = Dictionary{String, Int}()
    itd = Dictionary{String, Int}()
    ind = Dictionary{String, Int}()
    ied = Dictionary{String, Int}()
    for (e, n) in el
        if haskey(ELEMENTS[:ISOTOPES], e)
            get!(eld, e, 0)
            eld[e] += n
            get!(ied, e, 0)
            ied[e] += n
        elseif first(ELEMENTS[:ISOTOPES][get(ELEMENTS[:PARENTS], e, e)]) == e 
            p = get(ELEMENTS[:PARENTS], e, e)
            get!(ind, p, 0)
            ind[p] += n
            get!(eld, p, 0)
            eld[p] += n
            get!(ied, p, 0)
            ied[p] += n
        else 
            p = get(ELEMENTS[:PARENTS], e, e)
            get!(ind, e, 0)
            ind[e] += n
            get!(itd, e, 0)
            itd[e] += n
            get!(ied, e, 0)
            ied[e] += n
        end
    end
    # element => isoptope pairs
    # remove first
    ei = mapreduce(vcat, collect(keys(eld))) do e
        v = map(get(ELEMENTS[:ISOTOPES], e, e)) do x
            (e, x)
        end
        deleteat!(v, 1)
    end
    sort!(ei; by = x -> ELEMENTS[:ABUNDANCE][last(x)], rev = true)
    ab = isotopicabundance(el; ignore_isotopes = true)
    threshn = @match abtype begin
        :input => thresh * ab
        :max  => thresh * ab
        _     => thresh
    end
    # serve abundance as sums
    ite, ita = rec_isotopes!([ied], [abundance * ab], eld, itd, ind, ei, 1, abundance * ab, threshn)
    # Normalize after resolution check?
    ita = abtype == :max ? ita ./ maximum(ita) .* abundance : abtype == :input ? ita ./ first(ita) .* abundance : ita
    itm = map(mmi, ite, repeat([net_charge], length(ite)))
    id = sortperm(itm)
    filter!(x -> >=(ita[x], thresh), id)
    ite = ite[id]
    ita = ita[id]
    itm = itm[id]
    table ? Table(; Element = ite, Mass = itm, Abundance = ita) : ite
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
    # use [12C] for fix 
    for (e, n) in elements
        haskey(ELEMENTS[:PARENTS], e) || continue
        k = ELEMENTS[:PARENTS][e]
        haskey(ELEMENTS[:ISOTOPES], k) || continue
        get!(el, k, zeros(Int, length(ELEMENTS[:ISOTOPES][k])))
        haskey(ELEMENTS[:ISOTOPES], e) ? (el[e][begin] += n) : (ignore_isotopes || (el[k][findfirst(==(e), ELEMENTS[:ISOTOPES][k])] += n))
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

function initial_proportion(el::Dictionary, it::Dictionary, re::Dictionary, ir::Dictionary)
    p = 1
    ea = add_elements(el, re)
    ia = add_elements(it, ir)
    for e in keys(ea)
        n = 0
        for i in ISOTOPES[e]
            m = get(ia, i, 0)
            n += m
            p *= factorial(m)
        end
        p /= factorial(n + get(ea, e, 0), get(ea, e, 0))
    end
    for e in keys(el)
        n = 0
        for i in ISOTOPES[e]
            m = get(it, i, 0)
            n += m
            p /= factorial(m)
        end
        p *= factorial(n + get(el, e, 0), get(el, e, 0))
    end 
    for e in keys(re)
        n = 0
        for i in ISOTOPES[e]
            m = get(ir, i, 0)
            n += m
            p /= factorial(m)

        end
        p *= factorial(n + get(re, e, 0), get(re, e, 0))
    end 
    p
end

function update_proportion(prev_proportion, nold, nnew, delta)
    prev_proportion * factorial(nold, nold - delta) / factorial(nnew + delta, nnew)
end

function rec_isotopes!(ite::Vector, ita::Vector, el::Dictionary, it::Dictionary, ind::Dictionary, ei::Vector, ip::Int, prev_abundance, threshold)
    new_ip = -1
    for (e, i) in @views ei[ip:end]
        ne = get(el, e, 0)
        new_ip += 1
        if ne > get(ind, e, 0)
            new_el = deepcopy(el)
            new_it = deepcopy(it)
            new_el[e] -= 1
            get!(new_it, i, 0)
            ab = update_abundance(prev_abundance, e, i, el[e] - get(ind, e, 0), new_it[i] - get(ind, i, 0), 1)
            new_it[i] += 1
            ab < threshold && break
            push!(ite, add_elements(new_el, new_it))
            push!(ita, ab)
            rec_isotopes!(ite, ita, new_el, new_it, ind, ei, ip + new_ip, ab, threshold)
        end
    end
    ite, ita
end

# ab as proportion
function rec_2isotopes!(ite::Vector, ita::Vector, el::Dictionary, it::Dictionary, re::Dictionary, ir::Dictionary, ei::Vector, ip::Int, prev_proportion)
    new_ip = -1
    for (e, i) in @views ei[ip:end]
        ne = get(el, e, 0)
        nr = get(ir, i, 0)
        new_ip += 1
        if ne > 0 && nr > 0
            new_el = deepcopy(el)
            new_it = deepcopy(it)
            new_re = deepcopy(re)
            new_ir = deepcopy(ir)
            new_el[e] -= 1
            new_ir[i] -= 1
            get!(new_it, i, 0)
            get!(new_re, e, 0)
            ab = update_proportion(prev_proportion, el[e], new_it[i], 1)
            ab = update_proportion(ab, ir[i], new_re[e], 1)
            new_it[i] += 1
            new_re[e] += 1
            push!(ite, add_elements(new_el, new_it))
            push!(ita, ab)
            rec_2isotopes!(ite, ita, new_el, new_it, new_re, new_ir, ei, ip + new_ip, ab)
        end
    end
    ite, ita
end

function distribute_isotopes!(ep::Dictionary, el::Dictionary, it::Dictionary, re::Dictionary, ir::Dictionary)
    c = unique(map(e -> get(ELEMENTS[:PARENTS], e, e), collect(keys(ep))))
    els = [el]
    its = [it]
    res = [re]
    irs = [ir]
    for e in c 
        i = filter(i -> i != e && get(ep, i, 0) > 0, ELEMENTS[:ISOTOPES][e])
        nis = map(x -> ep[x], i)
        ni = sum(nis)
        if get(re, e, 0) >= ni 
            for re in res
                re[e] -= ni 
            end
            for ii in i 
                for ir in irs
                    get!(ir, ii, 0)
                    ir[ii] += ep[ii]
                end
            end
            continue
        end
        new_els = empty(els)
        new_its = empty(its)
        new_res = empty(res)
        new_irs = empty(irs)
        for ex in multiexponents(length(i), get(re, e, 0)) 
            any(x -> <(x...), zip(nis, ex)) && continue 
            delta = nis .- ex
            for ir in irs
                push!(new_irs, deepcopy(ir))
                for (x, ii) in zip(ex, i)
                    get!(last(new_irs), ii, 0)
                    last(new_irs)[ii] += x
                end
            end
            for it in its
                push!(new_its, deepcopy(it))
                for (d, ii) in zip(delta, i)
                    get!(last(new_its), ii, 0)
                    last(new_its)[ii] += d
                end
            end
        end
        for el in els
            for _ in eachindex(new_its)
                push!(new_els, deepcopy(el))
                get!(last(new_els), e, 0)
                last(new_els)[e] -= ni - get(re, e, 0)
            end
        end
        for re in res
            for _ in eachindex(new_its)
                push!(new_res, deepcopy(re))
                last(new_res)[e] = 0
            end
        end
        els = new_els
        its = new_its 
        res = new_res
        irs = new_irs
    end
    els, its, res, irs
end

# Distribute isotopes: setdiff without deleting overlapped/calculate 1 fragment * n
function isotopes_proportion(ep::Dictionary, ef::Dictionary, er::Dictionary, net_charge)
    eld = Dictionary{String, Int}()
    itd = Dictionary{String, Int}()
    for (e, n) in pairs(ef)
        if haskey(ELEMENTS[:ISOTOPES], e)
            get!(eld, e, 0)
            eld[e] += n
        else
            get!(itd, e, 0)
            itd[e] += n
        end
    end
    red = Dictionary{String, Int}()
    ird = Dictionary{String, Int}()
    for (e, n) in pairs(er)
        if haskey(ELEMENTS[:ISOTOPES], e)
            get!(red, e, 0)
            red[e] += n
        else
            throw(ArgumentError("MS2 gragment cannot contain isotopes."))
        end
    end
    # give isotope to red/ird first, multiple way ex 2C => C, 1 13C, 1 14C, can be C+13C => 14C or C+14C => 13C
    # for (e, n) in pairs(ep) 
    #     if get(ELEMENTS[:PARENTS], e, e) == e
    #         continue
    #     else
    #         p = get(ELEMENTS[:PARENTS], e, e)
    #         delta1 = min(n, get(red, p, 0))
    #         delta2 = n - delta1
    #         if delta1 > 0 
    #             red[p] -= delta1
    #             get!(ird, e, 0)
    #             ird[e] += delta1
    #         end 
    #         if delta2 > 0
    #             eld[p] -= delta2 
    #             get!(itd, e, 0)
    #             itd[e] += delta2 
    #         end
    #     end
    # end
    els, its, res, irs = distribute_isotopes!(ep, eld, itd, red, ird)
    # element => isoptope pairs
    # remove first
    ei = mapreduce(vcat, collect(keys(eld))) do e
        v = map(get(ELEMENTS[:ISOTOPES], e, e)) do x
            (e, x)
        end
        deleteat!(v, 1)
    end
    sort!(ei; by = x -> ELEMENTS[:ABUNDANCE][last(x)], rev = true)
    # initial proportion
    ite = []
    ita = []
    for (eld, itd, red, ird) in zip(els, its, res, irs)
        p = initial_proportion(eld, itd, red, ird) 
        itee, itaa = rec_2isotopes!([add_elements(eld, itd)], [p], eld, itd, red, ird, ei, 1, p)
        for (itx, ity) in zip(itee, itaa) 
            any(==(itx), ite) && continue
            push!(ite, itx)
            push!(ita, ity)
        end
    end
    itm = map(mmi, ite, repeat([net_charge], length(ite)))
    id = sortperm(itm)
    ite = ite[id]
    ita = ita[id]
    itm = itm[id]
    ite, itm, ita
end

function isobaric_sum(element::T, mmm, abundance, m_tol) where {T}
    e = T[]
    m = Float64[0]
    a = Vector{Float64}[]
    for (ee, mm, aa) in zip(element, mmm, abundance)
        if any(r -> in(mm, r), makecrit_delta(m_tol, last(m)))
            m[end] = (sum(last(a)) * m[end] + aa * mm) / (sum(last(a)) + aa)
            push!(last(e), ee)
            push!(last(a), aa)
        else
            push!(e, [ee])
            push!(a, [aa])
            push!(m, mm)
        end
    end
    popfirst!(m)
    for (x, y) in zip(e, a)
        ord = sortperm(y; rev = true)
        x .= x[ord]
        y .= y[ord]
    end
    e, m, a
end