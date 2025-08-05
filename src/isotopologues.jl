# Use Isotopes wrapper?
"""
    isotopologues(chemical::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6))
    isotopologues(formula::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6), net_charge = 0)
    isotopologues(chemicalpair::ChemicalPair, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6))) 
    isotopologues(formulapair::Pair{<: AbstractString, <: AbstractString}, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), net_charge = (1, 1))

Isotopologues of a single `chemical`, `formula`, or MS/MS precursor-product pairs (`chemicalpair`, and `formulapair`). Only isotopic abundance of parent elements are considered, and isotopes are viewed as intentionally labeled elements. 

* `abundance` sets the abundance of the isotope specified by `abtype`. 
* `abtype`
    * `:max`: the most abundant isotopologue.
    * `:input`: the input isotopologue.
    * `:total`: sum of total isotopologues.
* `threshold` can be a number or criteria, representing the lower limit of abundance. 
* `isobaric` determines whether groups isobars based on given tolerance. If the input is a subtype of `AbstractChemical`, it creates an `Isobars`; otherwise, it groups formulas as a vector. 
* `mz_tol` and `mm_tol` are tolerances of m/z and molecular mass for isobars.
    * number: acceptable range of absolute value, i.e. 0.01 for ±0.01
    * criteria: acceptable range of absolute and relative value, i.e. `crit(0.01, 20e-6)` for ±0.01 or 20 ppm error.
* `net_charge`: charges (positive or negative) of `formula`.

For MS/MS precursor-product pairs, `mz_tol`, `mm_tol`, and `net_charge` are pairs representing values for precursor and product, repectively. Product can also be neutral loss or chemical loss (`ChemicalLoss` or formula starting with `-`).

!!! Special precaution for applying to MS/MS precursor-product pairs
    Product must come from a single part or mutiple non-overlapping parts of precursor. Isobaric or isomeric products are not considered. For instance, 
    * PC 18:0/18:0 and fatty acid 18:0 is valid because two fatty acids are independent and identical. 
    * PC 18:0[D5]/18:0 and fatty acid 18:0[D5] is valid but the contribution of another fatty acid 18:0 is not considered and addional computation of this pair and summation with knowledge of fragmentation efficiency are required for the correct abundances. 
    * PC 18:1/18:0 and fatty acid 18:0 is also valid but requires additional computation of isobaric contribution of another fatty acid 18:1. 

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
    isotopologues_table(formulapair::Pair{<: AbstractString, <: AbstractString}, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), net_charge = (1, 1)) 
    isotopologues_table(tbl::Table; colchemical = :Chemical, colabundance = :Abundance, abundance = 1, colpreserve = setdiff(propertynames(tbl), [colchemical, colabundance]), kwargs...)
    isotopologues_table(v::Vector, abundance = 1; kwargs...)

A `Table` of isotopologues of a single `chemical`, `formula`, MS/MS precursor-product pairs (`chemicalpair`, and `formulapair`), or multiple chemicals in `v` or column `colchemical` of `tbl`. Only isotopic abundance of parent elements are considered, and isotopes are viewed as intentionally labeled elements. 

* `abundance` sets the abundance of the isotope specified by `abtype`. 
* `abtype`
    * `:max`: the most abundant isotopologue.
    * `:input`: the input isotopologue.
    * `:total`: sum of total isotopologues.
* `threshold` can be a number or criteria, representing the lower limit of abundance. 
* `isobaric` determines whether groups isobars based on given tolerance. If the input is a subtype of `AbstractChemical`, it creates an `Isobars`; otherwise, it groups formulas as a vector. 
* `mz_tol` and `mm_tol` are tolerances of m/z and molecular mass for isobars.
    * number: acceptable range of absolute value, i.e. 0.01 for ±0.01
    * criteria: acceptable range of absolute and relative value, i.e. `crit(0.01, 20e-6)` for ±0.01 or 20 ppm error.
* `net_charge`: charges (positive or negative) of `formula`.

For MS/MS precursor-product pairs, `mz_tol`, `mm_tol`, and `net_charge` are pairs representing values for precursor and product, repectively. Product can also be neutral loss or chemical loss (`ChemicalLoss` or formula starting with `-`).

!!! Special precaution for applying to MS/MS precursor-product pairs
    Product must come from a single part or mutiple non-overlapping parts of precursor. Isobaric or isomeric products are not considered. For instance, 
    * PC 18:0/18:0 and fatty acid 18:0 is valid because two fatty acids are independent and identical. 
    * PC 18:0[D5]/18:0 and fatty acid 18:0[D5] is valid but the contribution of another fatty acid 18:0 is not considered and addional computation of this pair and summation with knowledge of fragmentation efficiency are required for the correct abundances. 
    * PC 18:1/18:0 and fatty acid 18:0 is also valid but requires additional computation of isobaric contribution of another fatty acid 18:1. 

"""
isotopologues_table(cc::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6)) = 
    _isotopologues(cc, abundance; abtype, threshold, isobaric, mz_tol, mm_tol, table = true)

isotopologues_table(formula::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6), net_charge = 0) = 
    _isotopologues(formula, abundance; abtype, threshold, isobaric, mz_tol, mm_tol, net_charge, table = true)
    
isotopologues_table(cc::ChemicalPair, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6))) = 
    _isotopologues(cc, abundance; abtype, threshold, isobaric, mz_tol, mm_tol, table = true)

isotopologues_table(formula::Pair{<: AbstractString, <: AbstractString}, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), net_charge = (1, 1)) = 
    _isotopologues(formula, abundance; abtype, threshold, isobaric, mz_tol, mm_tol, net_charge, table = true)

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

Parent elements are viewed as major isotopes, and isotopic abundances of all elements are considered in computation. 
To compute isotopic abundance of chemicals with all isotopes labeled intentionally and not following natural distribution, set keyword argument `ignore_isotopes` true, and only parent elements are considered.  
"""
isotopicabundance(cc::AbstractChemical; ignore_isotopes = false) = isotopicabundance(chemicalformula(cc); ignore_isotopes)
isotopicabundance(formula::AbstractString; ignore_isotopes = false) = isotopicabundance(chemicalelements(formula); ignore_isotopes)
isotopicabundance(elements::Dictionary; ignore_isotopes = false) = isotopicabundance(pairs(elements); ignore_isotopes)
function isotopicabundance(elements::Union{<: Vector{<: Pair}, <: Dictionaries.PairDictionary}; ignore_isotopes = false)
    element_dict = Dictionary{String, Vector{Int}}()
    # use [12C] for fix 
    for (e, n) in elements
        haskey(ELEMENTS[:PARENTS], e) || continue
        k = ELEMENTS[:PARENTS][e]
        haskey(ELEMENTS[:ISOTOPES], k) || continue
        get!(element_dict, k, zeros(Int, length(ELEMENTS[:ISOTOPES][k])))
        haskey(ELEMENTS[:ISOTOPES], e) ? (element_dict[e][begin] += n) : (ignore_isotopes || (element_dict[k][findfirst(==(e), ELEMENTS[:ISOTOPES][k])] += n))
    end # use two vector?
    abundance_sum = 1
    for (e, ns) in pairs(element_dict)
        abundance = get.(Ref(ELEMENTS[:ABUNDANCE]), ELEMENTS[:ISOTOPES][e], 1)
        any(==(1), abundance) && continue
        n = sum(ns)
        id = sortperm(ns; rev = true)
        ns = ns[id]
        abundance = abundance[id]
        nspop1 = popfirst!(ns)
        f = try 
            factorial(n, nspop1)
        catch 
            factorial(big(n), nspop1)
        end
        abundance_sum *= f * popfirst!(abundance) ^ nspop1
        isempty(abundance) && continue
        abundance_sum *= mapreduce(*, ns, abundance) do x, y
            y ^ x / factorial(x)
        end
    end
    Float64(abundance_sum)
end

# ==========================================================================================================================
# Internal
function _isotopologues(formula::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6), net_charge = 0, table = true)
    if !isobaric
        it = _isotope_abundance(formula, abundance; abtype, threshold, net_charge)
        mass = it.Mass ./ max(1, abs(net_charge))
        chemical = map(chemicalformula, it.Element)
        table || return chemical
        if net_charge == 0 
            tbl = Table(; Isotopologues = chemical, Mass = mass, Abundance = it.Abundance) 
        else
            tbl = Table(; Isotopologues = chemical, MZ = mass, Abundance = it.Abundance) 
        end
        return tbl
    end
    it = _isotope_abundance(formula, abundance; abtype, threshold = threshold / 10, net_charge)
    # Table(; Formula = map(chemicalformula, it.Element), Mass = it.Mass, Abundance = it.Abundance)   
    m_tol = crit(net_charge == 0 ? mm_tol : mz_tol)
    it.Mass ./= max(1, abs(net_charge))
    element_single, mass_sum, abundance_single = isobaric_sum(it.Element, it.Mass, it.Abundance, m_tol)
    abundance_sum = map(sum, abundance_single)
    id = findall(>=(maximum(makecrit_value(crit(threshold), abundance))), abundance_sum)
    if table 
        if net_charge == 0 
            Table(; Isotopologues = map(x -> chemicalformula.(x), element_single[id]), Mass = mass_sum[id], Abundance = abundance_sum[id]) 
        else
            Table(; Isotopologues = map(x -> chemicalformula.(x), element_single[id]), MZ = mass_sum[id], Abundance = abundance_sum[id]) 
        end
    else 
        map(x -> chemicalformula.(x), element_single[id])
    end 
end

function _isotopologues(input_chemical::AbstractChemical, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = crit(0.01, 20e-6), mm_tol = crit(0.01, 20e-6), table = true)
    # if charge(input_chemical) == 0 
    #     it = _isotope_abundance(chemicalformula(input_chemical), abundance; abtype, threshold)
    #     return table ? Table(; Chemical = Isotopomers.(x, it.Element), Mass = it.Mass, Abundance = it.Abundance) : Isotopomers.(x, it.Element)
    # end
    # thresh / 10, to make error of abundance less than 10%
    if !isobaric
        it = _isotope_abundance(chemicalelements(input_chemical), abundance; abtype, threshold, net_charge = charge(input_chemical))
        mass = it.Mass ./ max(1, ncharge(input_chemical))
        chemical = Isotopomers.(input_chemical, it.Element)
        table || return chemical
        if charge(input_chemical) == 0 
            tbl = Table(; Isotopologues = chemical, Mass = mass, Abundance = it.Abundance) 
        else
            tbl = Table(; Isotopologues = chemical, MZ = mass, Abundance = it.Abundance) 
        end
        return tbl
    end
    it = _isotope_abundance(chemicalelements(input_chemical), abundance; abtype, threshold = threshold / 10, net_charge = charge(input_chemical))
    # it = Table(it; Element = map(x -> loss_elements!(unique_elements(x), adductelements(input_chemical)), it.Element))
    m_tol = crit(charge(input_chemical) == 0 ? mm_tol : mz_tol)
    it.Mass ./= max(1, ncharge(input_chemical))
    element_single, mass_sum, abundance_single = isobaric_sum(it.Element, it.Mass, it.Abundance, m_tol)
    abundance_sum = map(sum, abundance_single)
    id = findall(>=(maximum(makecrit_value(crit(threshold), abundance))), abundance_sum)
    if table 
        if charge(input_chemical) == 0 
            Table(; Isotopologues = map((x, y) -> Isobars(Isotopomers.(input_chemical, x), y), element_single[id], abundance_single[id]), Mass = mass_sum[id], Abundance = abundance_sum[id]) 
        else
            Table(; Isotopologues = map((x, y) -> Isobars(Isotopomers.(input_chemical, x), y), element_single[id], abundance_single[id]), MZ = mass_sum[id], Abundance = abundance_sum[id]) 
        end
    else 
        map((x, y) -> Isobars(Isotopomers.(input_chemical, x), y), element_single[id], abundance_single[id])
    end
end

function _isotopologues(formula::Pair{<: AbstractString, <: AbstractString}, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), net_charge = (1, 1), table = true)
    elements_precursor = chemicalelements(first(formula))
    element_precursor = filter(x -> haskey(ELEMENTS[:ISOTOPES], first(x)) && last(x) != 0, elements_precursor)
    isotope_precursor = filter(x -> !haskey(ELEMENTS[:ISOTOPES], first(x)) && last(x) != 0, elements_precursor)
    # any(x -> !haskey(ELEMENTS[:ISOTOPES], first(x)) && last(x) > 0, element_precursor) && throw(ArgumentError("Isotopologues table for MS/MS fragments of isotopic-labeled chemicals is not suppoeted"))
    if startswith(last(formula), "-")
        net_charge = (first(net_charge), first(net_charge) - last(net_charge))
        elements_product = chemicalelements(last(formula)[begin + 1:end])
        element_residual = unique_elements(filter(x -> haskey(ELEMENTS[:ISOTOPES], first(x)) && last(x) != 0, elements_product))
        isotope_residual = unique_elements(filter(x -> !haskey(ELEMENTS[:ISOTOPES], first(x)) && last(x) != 0, elements_product))
        element_product = loss_elements(element_precursor, element_residual)
        isotope_product = loss_elements(isotope_precursor, isotope_residual)
        loss = true
    else
        elements_product = chemicalelements(last(formula))
        element_product = unique_elements(filter(x -> haskey(ELEMENTS[:ISOTOPES], first(x)) && last(x) != 0, elements_product))
        isotope_product = unique_elements(filter(x -> !haskey(ELEMENTS[:ISOTOPES], first(x)) && last(x) != 0, elements_product))
        element_residual = loss_elements(element_precursor, element_product)
        loss = false
    end
    any(==(0), net_charge) && any(!=(0), net_charge) && throw(ArgumentError("Charges must be all non-zero or all zero."))
    # element_precursor_all, mass_precursor_all, abundance_precursor_all = _isotopologues(first(formula), abundance; abtype, threshold, isobaric, mz_tol = first(mz_tol), mm_tol = first(mm_tol), net_charge = first(net_charge), colision, table)
    if !isobaric
        it1 = _isotope_abundance(first(formula), abundance; abtype, threshold, net_charge = first(net_charge), normalize = false)
        it = _isotope_abundance_ms2(it1, isotope_precursor, element_product, isotope_product, element_residual, abundance; abtype, threshold, net_charge = last(net_charge))
        mass1 = it.Mass1 ./ max(1, abs(first(net_charge)))
        mass2 = it.Mass2 ./ max(1, abs(last(net_charge)))
        chemical = map(x -> chemicalformula(first(x)) => chemicalformula(last(x)), it.Element)
        table || return chemical
        if last(net_charge) == 0 
            tbl = Table(; Isotopologues = chemical, Mass1 = mass1, Mass2 = mass2, Abundance = it.Abundance) 
        else
            tbl = Table(; Isotopologues = chemical, MZ1 = mass1, MZ2 = mass2, Abundance = it.Abundance) 
        end
        return tbl
    end
    it1 = _isotope_abundance(first(formula), abundance; abtype, threshold = threshold / 10, net_charge = first(net_charge), normalize = false)
    it = _isotope_abundance_ms2(it1, isotope_precursor, element_product, isotope_product, element_residual, abundance; abtype, threshold, net_charge = last(net_charge))
    it.Mass1 ./= max(1, first(net_charge))
    it.Mass2 ./= max(1, last(net_charge))
    m_tol = (crit(first(net_charge) == 0 ? first(mm_tol) : first(mz_tol)), crit(last(net_charge) == 0 ? last(mm_tol) : last(mz_tol)))
    element_pair_single, mass_precursor_sum, mass_product_sum, abundance_pair_single = isobaric_sum(it.Element, it.Mass1, it.Mass2, it.Abundance, m_tol)
    abundance_pair_sum = map(sum, abundance_pair_single)
    id = findall(>=(maximum(makecrit_value(crit(threshold), abundance))), abundance_pair_sum)
    if table 
        if last(net_charge) == 0 
            Table(; Isotopologues = map(x -> map(y -> Pair(chemicalformula(first(y)), chemicalformula(last(y))), x), element_pair_single[id]), Mass1 = mass_precursor_sum[id], Mass2 = mass_product_sum[id], Abundance = abundance_pair_sum[id]) 
        else
            Table(; Isotopologues = map(x -> map(y -> Pair(chemicalformula(first(y)), chemicalformula(last(y))), x), element_pair_single[id]), MZ1 = mass_precursor_sum[id], MZ2 = mass_product_sum[id], Abundance = abundance_pair_sum[id]) 
        end
    else 
        map(x -> map(y -> Pair(chemicalformula(first(y)), chemicalformula(last(y))), x), element_pair_single[id])
    end 
end

function _isotopologues(cp::ChemicalPair, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), isobaric = true, mz_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), mm_tol = (crit(0.35, 5e-4), crit(0.01, 20e-6)), table = true)
    net_charge = (ncharge(cp.precursor), ncharge(cp.product))
    elements_precursor = chemicalelements(cp.precursor)
    elements_product = chemicalelements(cp.product)
    element_precursor = filter(x -> haskey(ELEMENTS[:ISOTOPES], first(x)) && last(x) != 0, elements_precursor)
    isotope_precursor = filter(x -> !haskey(ELEMENTS[:ISOTOPES], first(x)) && last(x) != 0, elements_precursor)
    # any(x -> !haskey(ELEMENTS[:ISOTOPES], first(x)) && last(x) > 0, element_precursor) && throw(ArgumentError("Isotopologues table for MS/MS fragments of isotopic-labeled chemicals is not suppoeted"))  
    if cp.product isa ChemicalLoss
        net_charge = (first(net_charge), first(net_charge) - last(net_charge))
        element_residual = unique_elements(filter(x -> haskey(ELEMENTS[:ISOTOPES], first(x)) && last(x) != 0, elements_product))
        isotope_residual = unique_elements(filter(x -> !haskey(ELEMENTS[:ISOTOPES], first(x)) && last(x) != 0, elements_product))
        element_product = loss_elements(element_precursor, element_residual)
        isotope_product = loss_elements(isotope_precursor, isotope_residual)
        loss = true
    else
        element_product = unique_elements(filter(x -> haskey(ELEMENTS[:ISOTOPES], first(x)) && last(x) != 0, elements_product))
        isotope_product = unique_elements(filter(x -> !haskey(ELEMENTS[:ISOTOPES], first(x)) && last(x) != 0, elements_product))
        element_residual = loss_elements(element_precursor, element_product)
        loss = false
    end
    any(==(0), net_charge) && any(!=(0), net_charge) && throw(ArgumentError("Charges must be all non-zero or all zero."))
    if !isobaric
        it1 = _isotope_abundance(chemicalelements(cp.precursor), abundance; abtype, threshold, net_charge = first(net_charge), normalize = false)
        it = _isotope_abundance_ms2(it1, isotope_precursor, element_product, isotope_product, element_residual, abundance; abtype, threshold, net_charge = last(net_charge))
        mass1 = it.Mass1 ./ max(1, abs(first(net_charge)))
        mass2 = it.Mass2 ./ max(1, abs(last(net_charge)))
        chemical = loss ? ChemicalPair.(Isotopomers.(cp.precursor, first.(it.Element)), Isotopomers.(cp.product, loss_elements.(first.(it.Element), last.(it.Element)))) : 
            ChemicalPair.(Isotopomers.(cp.precursor, first.(it.Element)), Isotopomers.(cp.product, last.(it.Element)))
        table || return chemical
        if last(net_charge) == 0 
            tbl = Table(; Isotopologues = chemical, Mass1 = mass1, Mass2 = mass2, Abundance = it.Abundance) 
        else
            tbl = Table(; Isotopologues = chemical, MZ1 = mass1, MZ2 = mass2, Abundance = it.Abundance) 
        end
        return tbl
    end
    it1 = _isotope_abundance(chemicalelements(cp.precursor), abundance; abtype, threshold = threshold / 10, net_charge = first(net_charge), normalize = false)
    it = _isotope_abundance_ms2(it1, isotope_precursor, element_product, isotope_product, element_residual, abundance; abtype, threshold, net_charge = last(net_charge))
    it.Mass1 ./= max(1, first(net_charge))
    it.Mass2 ./= max(1, last(net_charge))
    m_tol = (crit(first(net_charge) == 0 ? first(mm_tol) : first(mz_tol)), crit(last(net_charge) == 0 ? last(mm_tol) : last(mz_tol)))
    element_pair_single, mass_precursor_sum, mass_product_sum, abundance_pair_single = isobaric_sum(it.Element, it.Mass1, it.Mass2, it.Abundance, m_tol)
    abundance_pair_sum = map(sum, abundance_pair_single)
    id = findall(>=(maximum(makecrit_value(crit(threshold), abundance))), abundance_pair_sum)
    if table 
        if last(net_charge) == 0 && loss
            Table(; Isotopologues = map((x, a) -> Isobars(map(y -> ChemicalPair(Isotopomers(cp.precursor, first(y)), Isotopomers(cp.product, loss_elements(first(y), last(y)))), x), a), element_pair_single[id], abundance_pair_single[id]), Mass1 = mass_precursor_sum[id], Mass2 = mass_product_sum[id], Abundance = abundance_pair_sum[id]) 
        elseif last(net_charge) == 0
            Table(; Isotopologues = map((x, a) -> Isobars(map(y -> ChemicalPair(Isotopomers(cp.precursor, first(y)), Isotopomers(cp.product, last(y))), x), a), element_pair_single[id], abundance_pair_single[id]), Mass1 = mass_precursor_sum[id], Mass2 = mass_product_sum[id], Abundance = abundance_pair_sum[id]) 
        elseif loss
            Table(; Isotopologues = map((x, a) -> Isobars(map(y -> ChemicalPair(Isotopomers(cp.precursor, first(y)), Isotopomers(cp.product, loss_elements(first(y), last(y)))), x), a), element_pair_single[id], abundance_pair_single[id]), MZ1 = mass_precursor_sum[id], MZ2 = mass_product_sum[id], Abundance = abundance_pair_sum[id]) 
        else
            Table(; Isotopologues = map((x, a) -> Isobars(map(y -> ChemicalPair(Isotopomers(cp.precursor, first(y)), Isotopomers(cp.product, last(y))), x), a), element_pair_single[id], abundance_pair_single[id]), MZ1 = mass_precursor_sum[id], MZ2 = mass_product_sum[id], Abundance = abundance_pair_sum[id]) 
        end
    elseif loss
        map((x, a) -> Isobars(map(y -> ChemicalPair(Isotopomers(cp.precursor, first(y)), Isotopomers(cp.product, loss_elements(first(y), last(y)))), x), a), element_pair_single[id], abundance_pair_single[id])
    else
        map((x, a) -> Isobars(map(y -> ChemicalPair(Isotopomers(cp.precursor, first(y)), Isotopomers(cp.product, last(y))), x), a), element_pair_single[id], abundance_pair_single[id])
    end 
end

# loewest level
_isotope_abundance(x::AbstractString, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), net_charge = 0, table = true, normalize = true) = 
    _isotope_abundance(chemicalelements(x), abundance; abtype, threshold, net_charge, table, normalize)
function _isotope_abundance(input_element::Vector, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), net_charge = 0, table = true, normalize = true)
    abundance_cutoff = maximum(makecrit_value(crit(threshold), abundance))
    # record elements change
    element_dictionary = Dictionary{String, Int}()
    isotope_dictionary = Dictionary{String, Int}()
    fix_dictionary = Dictionary{String, Int}()
    first_element_dictionary = Dictionary{String, Int}()
    for (e, n) in input_element
        if haskey(ELEMENTS[:ISOTOPES], e)
            get!(element_dictionary, e, 0)
            element_dictionary[e] += n
            get!(first_element_dictionary, e, 0)
            first_element_dictionary[e] += n
        elseif first(ELEMENTS[:ISOTOPES][get(ELEMENTS[:PARENTS], e, e)]) == e 
            p = get(ELEMENTS[:PARENTS], e, e)
            get!(fix_dictionary, p, 0)
            fix_dictionary[p] += n
            get!(element_dictionary, p, 0)
            element_dictionary[p] += n
            get!(first_element_dictionary, p, 0)
            first_element_dictionary[p] += n
        else 
            p = get(ELEMENTS[:PARENTS], e, e)
            get!(fix_dictionary, e, 0)
            fix_dictionary[e] += n
            get!(isotope_dictionary, e, 0)
            isotope_dictionary[e] += n
            get!(first_element_dictionary, e, 0)
            first_element_dictionary[e] += n
        end
    end
    # element => isoptope pairs
    # remove first
    element_isotope_pair = mapreduce(vcat, collect(keys(element_dictionary))) do e
        v = map(get(ELEMENTS[:ISOTOPES], e, e)) do x
            (e, x)
        end
        deleteat!(v, 1)
    end
    sort!(element_isotope_pair; by = x -> ELEMENTS[:ABUNDANCE][last(x)], rev = true)
    first_proportion = isotopicabundance(input_element; ignore_isotopes = true)
    abundance_cutoff_normalize = @match abtype begin
        :input => abundance_cutoff * first_proportion
        :max  => abundance_cutoff * first_proportion
        _     => abundance_cutoff
    end
    # serve abundance as sums
    element_chemical, abundance_chemical = rec_addisotopes!([first_element_dictionary], [abundance * first_proportion], element_dictionary, isotope_dictionary, isotope_dictionary, element_isotope_pair, 1, abundance * first_proportion, abundance_cutoff_normalize)
    # Normalize after resolution check?
    abundance_chemical_normalize = abtype == :max ? abundance_chemical ./ maximum(abundance_chemical) .* abundance : abtype == :input ? abundance_chemical ./ first(abundance_chemical) .* abundance : abundance_chemical
    mass_chemical = map(mmi, element_chemical, repeat([net_charge], length(element_chemical)))
    id = sortperm(mass_chemical)
    filter!(x -> >=(abundance_chemical_normalize[x], abundance_cutoff), id)
    element_chemical = element_chemical[id]
    abundance_chemical = normalize ? abundance_chemical_normalize[id] : abundance_chemical[id]
    mass_chemical = mass_chemical[id]
    table ? Table(; Element = element_chemical, Mass = mass_chemical, Abundance = abundance_chemical) : element_chemical
end

function _isotope_abundance_ms2(it1, isotope_precursor, element_product, isotope_product, element_residual, abundance = 1; abtype = :max, threshold = crit(abundance * 1e-4, 1e-4), net_charge = 0, normalize = true, table = true)
    abundance_cutoff = maximum(makecrit_value(crit(threshold), abundance))
    element_pair = Pair[]
    mass_product = Float64[]
    mass_precursor = Float64[]
    abundance_pair = Float64[]
    msfix = mmi(isotope_product, net_charge)
    for (element_precursor_i, mass_precursor_i, abundance_precursor_i) in zip(it1.Element, it1.Mass, it1.Abundance)
        element_product_is, mass_product_is, proportion_product_is = isotopes_proportion(loss_elements(element_precursor_i, isotope_precursor), element_product, element_residual, msfix, net_charge)
        id = sortperm(mass_product_is)
        element_pair = vcat(element_pair, Ref(element_precursor_i) .=> add_elements.(Ref(isotope_product), element_product_is[id]))
        mass_product = vcat(mass_product, mass_product_is[id])
        mass_precursor = vcat(mass_precursor, repeat([mass_precursor_i], length(element_product_is[id])))
        abundance_pair = vcat(abundance_pair, abundance_precursor_i .* proportion_product_is[id])
    end
    abundance_pair_normalize = abtype == :max ? abundance_pair ./ maximum(abundance_pair) .* abundance : abtype == :input ? abundance_pair ./ first(abundance_pair) .* abundance : abundance_pair
    id = findall(>=(abundance_cutoff), abundance_pair_normalize)
    element_pair = element_pair[id]
    mass_precursor = mass_precursor[id]
    mass_product = mass_product[id]
    abundance_pair = normalize ? abundance_pair_normalize[id] : abundance_pair[id]
    table ? Table(; Element = element_pair, Mass1 = mass_precursor, Mass2 = mass_product, Abundance = abundance_pair) : element_pair
end

function rec_addisotopes!(element_vec::Vector, abundance_vec::Vector, element_dictionary::Dictionary, isotope_dictionary::Dictionary, fix_dictionary::Dictionary, element_isotope_pair::Vector, isotope_position::Int, prev_abundance, threshold)
    new_isotope_position = -1
    for (e, i) in @views element_isotope_pair[isotope_position:end]
        ne = get(element_dictionary, e, 0)
        new_isotope_position += 1
        if ne > get(fix_dictionary, e, 0)
            new_element_dictionary = deepcopy(element_dictionary)
            new_isotope_dictionary = deepcopy(isotope_dictionary)
            new_element_dictionary[e] -= 1
            get!(new_isotope_dictionary, i, 0)
            abundance = update_abundance(prev_abundance, e, i, element_dictionary[e] - get(fix_dictionary, e, 0), new_isotope_dictionary[i] - get(fix_dictionary, i, 0), 1)
            new_isotope_dictionary[i] += 1
            abundance < threshold && break
            push!(element_vec, add_elements(new_element_dictionary, new_isotope_dictionary))
            push!(abundance_vec, abundance)
            rec_addisotopes!(element_vec, abundance_vec, new_element_dictionary, new_isotope_dictionary, fix_dictionary, element_isotope_pair, isotope_position + new_isotope_position, abundance, threshold)
        end
    end
    element_vec, abundance_vec
end

function update_abundance(prev_abundance, old_element, new_element, nold, nnew, delta)
    x = get(ELEMENTS[:ABUNDANCE], old_element, 1)
    y = get(ELEMENTS[:ABUNDANCE], new_element, 1)
    (x == 1 || y == 1) && return prev_abundance
    prev_abundance * factorial(nold, nold - delta) / factorial(nnew + delta, nnew) * (y / x) ^ delta
end

# Distribute isotopes: setdiff without deleting overlapped/calculate 1 fragment * n
function isotopes_proportion(element_precursor_dictionary::Dictionary, element_product_dictionary::Dictionary, element_residual_dictionary::Dictionary, msfix, net_charge)
    first_isotope_product_dictionary = Dictionary{String, Int}()
    first_isotope_residual_dictionary = Dictionary{String, Int}()
    first_element_product_dictionary = deepcopy(element_product_dictionary)
    first_element_residual_dictionary = deepcopy(element_residual_dictionary)
    # element => isoptope pairs
    # remove first
    element_isotope_pair = mapreduce(vcat, collect(keys(element_product_dictionary))) do e
        v = map(get(ELEMENTS[:ISOTOPES], e, e)) do x
            (e, x)
        end
        deleteat!(v, 1)
    end
    sort!(element_isotope_pair; by = x -> ELEMENTS[:ABUNDANCE][last(x)], rev = true)
    element_product_vec, isotope_product_vec, element_residual_vec, isotope_residual_vec = distribute_isotopes!(element_precursor_dictionary, first_element_product_dictionary, first_isotope_product_dictionary, first_element_residual_dictionary, first_isotope_residual_dictionary)
    # initial proportion
    element_pair = Dictionary[]
    proportion_pair = Float64[]
    for (element_product, isotope_product, element_residual, isotope_residual) in zip(element_product_vec, isotope_product_vec, element_residual_vec, isotope_residual_vec)
        proportion = initial_proportion(element_product, isotope_product, element_residual, isotope_residual)
        element_pair_vec, proportion_pair_vec = rec_moveisotopes!([add_elements(element_product, isotope_product)], [proportion], element_product, isotope_product, element_residual, isotope_residual, element_isotope_pair, 1, proportion)
        for (element_pairi, proportion_pairi) in zip(element_pair_vec, proportion_pair_vec) 
            any(==(element_pairi), element_pair) && continue
            push!(element_pair, element_pairi)
            push!(proportion_pair, proportion_pairi)
        end
    end
    mass_pair = map(mmi, element_pair, repeat([net_charge], length(element_pair))) .+ msfix
    id = sortperm(mass_pair)
    element_pair = element_pair[id]
    proportion_pair = proportion_pair[id]
    mass_pair = mass_pair[id]
    element_pair, mass_pair, proportion_pair
end

function distribute_isotopes!(element_precursor_dictionary::Dictionary, element_product_dictionary::Dictionary, isotope_product_dictionary::Dictionary, element_residual_dictionary::Dictionary, isotope_residual_dictionary::Dictionary)
    elements = unique(map(e -> get(ELEMENTS[:PARENTS], e, e), collect(keys(element_precursor_dictionary))))
    element_product_vec = [element_product_dictionary]
    isotope_product_vec = [isotope_product_dictionary]
    element_residual_vec = [element_residual_dictionary]
    isotope_residual_vec = [isotope_residual_dictionary]
    for e in elements
        isotopes = filter(i -> i != e && get(element_precursor_dictionary, i, 0) > 0, ELEMENTS[:ISOTOPES][e])
        nisotopes_vec = map(x -> element_precursor_dictionary[x], isotopes)
        nisotopes_total = sum(nisotopes_vec)
        if get(last(element_residual_vec), e, 0) >= nisotopes_total 
            for element_residual in element_residual_vec
                get!(element_residual, e, 0)
                element_residual[e] -= nisotopes_total 
            end
            for isotope in isotopes
                for isotope_residual in isotope_residual_vec
                    get!(isotope_residual, isotope, 0)
                    isotope_residual[isotope] += element_precursor_dictionary[isotope]
                end
            end
            continue
        end
        new_element_product_vec = empty(element_product_vec)
        new_isotope_product_vec = empty(isotope_product_vec)
        new_element_residual_vec = empty(element_residual_vec)
        new_isotope_residual_vec = empty(isotope_residual_vec)
        for nisotopes_residual_vec in multiexponents(length(isotopes), get(last(element_residual_vec), e, 0)) 
            any(x -> <(x...), zip(nisotopes_vec, nisotopes_residual_vec)) && continue 
            delta = nisotopes_vec .- nisotopes_residual_vec
            for isotope_residual in isotope_residual_vec
                push!(new_isotope_residual_vec, deepcopy(isotope_residual))
                for (x, isotope) in zip(nisotopes_residual_vec, isotopes)
                    get!(last(new_isotope_residual_vec), isotope, 0)
                    last(new_isotope_residual_vec)[isotope] += x
                end
            end
            for isotope_product in isotope_product_vec
                push!(new_isotope_product_vec, deepcopy(isotope_product))
                for (d, isotope) in zip(delta, isotopes)
                    get!(last(new_isotope_product_vec), isotope, 0)
                    last(new_isotope_product_vec)[isotope] += d
                end
            end
        end
        for element_product in element_product_vec
            for _ in eachindex(new_isotope_product_vec)
                push!(new_element_product_vec, deepcopy(element_product))
                get!(last(new_element_product_vec), e, 0)
                last(new_element_product_vec)[e] -= nisotopes_total - get(last(element_residual_vec), e, 0) - get(last(isotope_residual_vec), e, 0)
            end
        end
        for element_residual in element_residual_vec
            for _ in eachindex(new_isotope_product_vec)
                push!(new_element_residual_vec, deepcopy(element_residual))
                get!(last(new_element_residual_vec), e, 0)
                last(new_element_residual_vec)[e] = 0
            end
        end
        element_product_vec = new_element_product_vec
        isotope_product_vec = new_isotope_product_vec 
        element_residual_vec = new_element_residual_vec
        isotope_residual_vec = new_isotope_residual_vec
    end
    element_product_vec, isotope_product_vec, element_residual_vec, isotope_residual_vec
end

function initial_proportion(element_product::Dictionary, isotope_product::Dictionary, element_residual::Dictionary, isotope_residual::Dictionary)
    p = 1
    element_total = add_elements(element_product, element_residual)
    isotope_total = add_elements(isotope_product, isotope_residual)
    for e in keys(element_total)
        n = 0
        for i in ISOTOPES[e]
            m = get(isotope_total, i, 0)
            n += m
            p *= factorial(m)
        end
        p /= factorial(n + get(element_total, e, 0), get(element_total, e, 0))
    end
    for e in keys(element_product)
        n = 0
        for i in ISOTOPES[e]
            m = get(isotope_product, i, 0)
            n += m
            p /= factorial(m)
        end
        p *= factorial(n + get(element_product, e, 0), get(element_product, e, 0))
    end 
    for e in keys(element_residual)
        n = 0
        for i in ISOTOPES[e]
            m = get(isotope_residual, i, 0)
            n += m
            p /= factorial(m)

        end
        p *= factorial(n + get(element_residual, e, 0), get(element_residual, e, 0))
    end 
    p
end

# ab as proportion
function rec_moveisotopes!(element_vec::Vector, abundance_vec::Vector, element_product_dictionary::Dictionary, isotope_product_dictionary::Dictionary, element_residual_dictionary::Dictionary, isotope_residual_dictionary::Dictionary, element_isotope_pair::Vector, isotope_position::Int, prev_proportion)
    new_isotope_position = -1
    for (e, i) in @views element_isotope_pair[isotope_position:end]
        np = get(element_product_dictionary, e, 0)
        nr = get(isotope_residual_dictionary, i, 0)
        new_isotope_position += 1
        if np > 0 && nr > 0
            new_element_product_dictionary = deepcopy(element_product_dictionary)
            new_isotope_product_dictionary = deepcopy(isotope_product_dictionary)
            new_element_residual_dictionary = deepcopy(element_residual_dictionary)
            new_isotope_residual_dictionary = deepcopy(isotope_residual_dictionary)
            new_element_product_dictionary[e] -= 1
            new_isotope_residual_dictionary[i] -= 1
            get!(new_isotope_product_dictionary, i, 0)
            get!(new_element_residual_dictionary, e, 0)
            proportion = update_proportion(prev_proportion, element_product_dictionary[e], new_isotope_product_dictionary[i], 1)
            proportion = update_proportion(proportion, isotope_residual_dictionary[i], new_element_residual_dictionary[e], 1)
            new_isotope_product_dictionary[i] += 1
            new_element_residual_dictionary[e] += 1
            push!(element_vec, add_elements(new_element_product_dictionary, new_isotope_product_dictionary))
            push!(abundance_vec, proportion)
            rec_moveisotopes!(element_vec, abundance_vec, new_element_product_dictionary, new_isotope_product_dictionary, new_element_residual_dictionary, new_isotope_residual_dictionary, element_isotope_pair, isotope_position + new_isotope_position, proportion)
        end
    end
    element_vec, abundance_vec
end

function update_proportion(prev_proportion, nold, nnew, delta)
    prev_proportion * factorial(nold, nold - delta) / factorial(nnew + delta, nnew)
end

function isobaric_sum(element::T, mass::Vector{<: Number}, abundance, m_tol) where {T}
    element_single = T[]
    mass_sum = Float64[0]
    abundance_single = Vector{Float64}[]
    for (e, m, a) in zip(element, mass, abundance)
        if any(r -> in(m, r), makecrit_delta(m_tol, last(mass_sum)))
            mass_sum[end] = (sum(last(abundance_single)) * mass_sum[end] + a * m) / (sum(last(abundance_single)) + a)
            push!(last(element_single), e)
            push!(last(abundance_single), a)
        else
            push!(element_single, [e])
            push!(abundance_single, [a])
            push!(mass_sum, m)
        end
    end
    popfirst!(mass_sum)
    for (x, y) in zip(element_single, abundance_single)
        ord = sortperm(y; rev = true)
        x .= x[ord]
        y .= y[ord]
    end
    element_single, mass_sum, abundance_single
end

function isobaric_sum(element::T, mass1::Vector{<: Number}, mass2::Vector{<: Number}, abundance, m_tol) where {T}
    element_single = Vector{T}[]
    mass1_sum = Vector{Float64}[[0]]
    mass2_sum = Vector{Float64}[[0]]
    abundance_single = Vector{Vector{Float64}}[]
    for (e, m1, m2, a) in zip(element, mass1, mass2, abundance)
        push_m1 = false
        push_m2 = false
        for (i, (rm1, rm2)) in enumerate(zip(last(mass1_sum), last(mass2_sum)))
            if any(r -> in(m1, r), makecrit_delta(first(m_tol), rm1)) 
                push_m1 = true
                if any(r -> in(m2, r), makecrit_delta(first(m_tol), rm2))
                    mass1_sum[end][i] = (sum(last(abundance_single)[i]) * mass1_sum[end][i] + a * m1) / (sum(last(abundance_single)[i]) + a)
                    mass2_sum[end][i] = (sum(last(abundance_single)[i]) * mass2_sum[end][i] + a * m2) / (sum(last(abundance_single)[i]) + a)
                    push!(last(element_single)[i], e)
                    push!(last(abundance_single)[i], a)
                    push_m2 = true
                    break 
                end
            end
        end
        if push_m1 && !push_m2
            push!(last(mass1_sum), m1)
            push!(last(mass2_sum), m2)
            push!(last(element_single), [e])
            push!(last(abundance_single), [a])
            id = sortperm(last(mass2_sum))
            mass1_sum[end] = mass1_sum[end][id]
            mass2_sum[end] = mass2_sum[end][id]
            element_single[end] = element_single[end][id]
            abundance_single[end] = abundance_single[end][id]
        elseif !push_m1
            push!(mass1_sum, [m1])
            push!(mass2_sum, [m2])
            push!(element_single, [[e]])
            push!(abundance_single, [[a]])
        end
    end
    popfirst!(mass1_sum)
    popfirst!(mass2_sum)
    mass1_sum = vcat(mass1_sum...)
    mass2_sum = vcat(mass2_sum...)
    element_single = vcat(element_single...)
    abundance_single = vcat(abundance_single...)
    for (x, y) in zip(element_single, abundance_single)
        ord = sortperm(y; rev = true)
        x .= x[ord]
        y .= y[ord]
    end
    element_single, mass1_sum, mass2_sum, abundance_single
end
