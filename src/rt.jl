"""
    isobars_rt(ion::AbstractChemical, lib::AbstractVector{T}; fwhm_exp = 0.2, fwhm_lib = 0.2, mz_tol = crit(0.01, 20e-6))
    isobars_rt(ion::AbstractChemical, lib::Table; libchemical = :AdductIon, libmz = :MZ, librt = :RT, libfwhm = :FWHM, mz_tol = crit(0.01, 20e-6))

Return a subvector of `lib` such that all elements are isobaric to `ion` within specific tolerance of mz and FWHM values. 

`fwhm_exp` and `fwhm_lib` can be a number or a vector matching the length of `exp` and `lib`. 

Keyword arguments `libchemical`, `libmz`, `librt`, and `libfwhm` specify the columns of `AdductIon` objects, m/z values, rt values, and FWHM. 

Tolerance can be a number or criteria, representing maximum allowed difference. For criteria with multiple allowed intervals, lying in any of them is considered to fulfill the criteria. 
"""
function isobars_rt(ion::AbstractChemical, lib::AbstractVector{T}; fwhm_exp = 0.2, fwhm_lib = nothing, mz_tol = crit(0.01, 20e-6)) where {T <: AbstractChemical}
    isnan(rt(ion)) && return T[]
    rrt = rt(ion)
    rmz1 = mz(ion)
    mz_tol = union(makecrit_delta(crit(mz_tol), rmz1)...)
    rt_tol = makerttol(fwhm_exp, fwhm_lib, rrt, lib)
    c = T[]
    for (i, rt_tol_) in zip(lib, rt_tol)
        ischemicalequal(ion, i) && continue
        irt = rt(i)
        isnan(irt) && continue
        in(irt, rt_tol_) || continue
        imz1 = mz(i)
        in(imz1, mz_tol) || continue
        push!(c, i)
    end
    c
end

"""
    isobars_rt_table(ion::AbstractChemical, lib::AbstractVector{T}; fwhm_exp = 0.2, fwhm_lib = 0.2, mz_tol = crit(0.01, 20e-6))
    isobars_rt_table(exp::Vector, lib::AbstractVector{T}; kwargs...)
    isobars_rt_table(exp::Table, lib::AbstractVector{T}; fwhm_lib = 0.2, expid = true, expchemical = :AdductIon, expmz = :MZ, exprt = :RT, expfwhm = :FWHM, mz_tol = crit(0.01, 20e-6))
    isobars_rt_table(ion::AbstractChemical, lib::Table; fwhm_exp = 0.2, libid = true, libchemical = :AdductIon, libmz = :MZ, librt = :RT, libfwhm = :FWHM, mz_tol = crit(0.01, 20e-6))
    isobars_rt_table(exp::Vector, lib::Table; fwhm_exp = 0.2, libid = true, libchemical = :AdductIon, libmz = :MZ, librt = :RT, libfwhm = :FWHM, mz_tol = crit(0.01, 20e-6))
    isobars_rt_table(exp::Table, lib::Table; expid = true, expchemical = :AdductIon, expmz = :MZ, exprt = :RT, expfwhm = :FWHM, libid = true, libchemical = :AdductIon, libmz = :MZ, librt = :RT, libfwhm = :FWHM, mz_tol = crit(0.01, 20e-6))

Return a `Table` such that all rows representing an ion from `lib` which are isobaric to `ion` or ions from `exp` within specific tolerance of mz and FWHM values. 

Tolerance can be a number or criteria, representing maximum allowed difference. For criteria with multiple allowed intervals, lying in any of them is considered to fulfill the criteria. 

`fwhm_exp` and `fwhm_lib` can be a number or a vector matching the length of `exp` and `lib`. 

Keyword arguments `libchemical`, `libmz`, `librt`, `libfwhm` specify the columns of `AdductIon` objects, m/z values, rt values, and FWHM from `lib`. 

Keyword arguments `expchemical`, `expmz`, `exprt`, and `expfwhm` specify the columns of `AdductIon` objects, m/z values, rt values, and FWHM from `exp`. 

Keyword arguments `libid` and `expid` determine whether creates columns to store row id of `lib` or `exp`.
"""
isobars_rt_table(ion::AbstractChemical, lib::AbstractVector{T}; fwhm_exp = 0.2, fwhm_lib = nothing, mz_tol = crit(0.01, 20e-6)) where {T <: AbstractChemical} = 
    isobars_rt_table(ion, rt(ion), mz(ion), lib; fwhm_exp, fwhm_lib, mz_tol)
function isobars_rt_table(ion::AbstractChemical, rrt, rmz1, lib::AbstractVector{T}; fwhm_exp = 0.2, fwhm_lib = nothing, mz_tol = crit(0.01, 20e-6)) where {T <: AbstractChemical}
    isnan(rrt) && return nothing
    mz_tol = union(makecrit_delta(crit(mz_tol), rmz1)...)
    rt_tol = makerttol(fwhm_exp, fwhm_lib, rrt, lib)
    c = T[]
    Δrt = Float64[]
    Δrtp = Float64[]
    Δmz = Float64[]
    Δmzp = Float64[]
    for (i, rt_tol_) in zip(lib, rt_tol)
        ischemicalequal(ion, i) && continue
        irt = rt(i)
        isnan(irt) && continue
        in(irt, rt_tol_) || continue
        imz1 = mz(i)
        in(imz1, mz_tol) || continue
        push!(c, i)
        push!(Δrt, rrt - irt)
        push!(Δrtp, last(Δrt) / irt * 100)
        push!(Δmz, rmz1 - imz1)
        push!(Δmzp, last(Δmz) / imz1 * 1e6)
    end
    id = sortperm(Δrt; by = abs)
    c = c[id]
    Δrt = Δrt[id]
    Δrtp = Δrtp[id]
    Δmz = Δmz[id]
    Δmzp = Δmzp[id]
    Table(; Isobar = c, ΔRT = Δrt, var"ΔRT(%)" = Δrtp, ΔMZ = Δmz, var"ΔMZ(ppm)" = Δmzp)
end

function isobars_rt_table(exp::Vector, lib::AbstractVector{T}; fwhm_exp = 0.2, fwhm_lib = nothing, kwargs...) where {T <: AbstractChemical}
    fwhm_exp = collecfwhm(fwhm_exp, exp)
    if length(exp) < Threads.nthreads()
        mapreduce(vcat, exp, fwhm_exp) do x, y
            it = isobars_rt_table(x, lib; fwhm_exp = y, fwhm_lib, kwargs...)
            Table((AdductIon = vectorize(x, length(it)), ), it)
        end |> Table
    else
        t = Vector{Table}(undef, length(exp))
        Threads.@threads for i in eachindex(t)
            it = isobars_rt_table(exp[i], lib; fwhm_exp = fwhm_exp[i], fwhm_lib, kwargs...)
            t[i] = Table((AdductIon = vectorize(exp[i], length(it)), ), it)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end

function isobars_rt_table(exp::Table, lib::AbstractVector{T}; fwhm_lib = nothing, expid = true, expchemical = :AdductIon, expmz = :MZ, exprt = :RT, expfwhm = :FWHM, mz_tol = crit(0.01, 20e-6)) where {T <: AbstractChemical}
    expchemical = Symbol(expchemical)
    expmz = isnothing(expmz) ? nothing : Symbol(expmz)
    exprt = isnothing(exprt) ? nothing : Symbol(exprt)
    expfwhm = isnothing(expfwhm) ? nothing : Symbol(expfwhm)
    rrt = isnothing(exprt) ? rt.(getproperty(exp, expchemical)) : getproperty(exp, exprt)
    rmz1 = isnothing(expmz) ? mz.(getproperty(exp, expchemical)) : getproperty(exp, expmz)
    fwhm_exp = isnothing(expfwhm) ? [0.2 for _ in eachindex(exp)] : getproperty(exp, expfwhm)
    if length(exp) < Threads.nthreads()
        mapreduce(vcat, eachindex(exp), getproperty(exp, expchemical), fwhm_exp, rrt, rmz1) do i, x, y, xrt, xmz1
            it = isobars_rt_table(x, xrt, xmz1, lib; fwhm_exp = y, fwhm_lib, mz_tol)
            Table(expid ? (ExpID = vectorize(i, length(it)), AdductIon = vectorize(x, length(it)), ) : (AdductIon = vectorize(x, length(it)), ), it)
        end |> Table
    else
        t = Vector{Table}(undef, length(exp))
        Threads.@threads for i in eachindex(t)
            it = isobars_rt_table(getproperty(exp, expchemical)[i], rrt[i], rmz1[i], lib; fwhm_exp = fwhm_exp[i], fwhm_lib, mz_tol)
            t[i] = Table(expid ? (ExpID = vectorize(i, length(it)), AdductIon = vectorize(getproperty(exp, expchemical)[i], length(it)), ) : (AdductIon = vectorize(getproperty(exp, expchemical)[i], length(id)), ), it)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end

function isobars_rt(ion::AbstractChemical, lib::Table; fwhm_exp = 0.2, libchemical = :AdductIon, libmz = :MZ, librt = :RT, libfwhm = :FWHM, mz_tol = crit(0.01, 20e-6))
    libchemical = Symbol(libchemical)
    libmz = isnothing(libmz) ? nothing : Symbol(libmz)
    librt = isnothing(librt) ? nothing : Symbol(librt)
    libfwhm = isnothing(libfwhm) ? nothing : Symbol(libfwhm)
    isnan(rt(ion)) && return eltype(getproperty(lib, libchemical))[]
    getmz = isnothing(libmz) ? (i -> mz(getproperty(i, libchemical))) : (i -> getproperty(i, libmz))
    getrt = isnothing(librt) ? (i -> rt(getproperty(i, libchemical))) : (i -> getproperty(i, librt))
    fwhm = isnothing(libfwhm) ? [fwhm_exp for _ in eachindex(lib)] : (getproperty(lib, libfwhm) .+ fwhm_exp) ./ 2
    rrt = rt(ion)
    rmz1 = mz(ion)
    mz_tol = union(makecrit_delta(crit(mz_tol), rmz1)...)
    rt_tol = [union(makecrit_delta(crit(x), rrt)...) for x in fwhm]
    c = eltype(getproperty(lib, libchemical))[]
    for (i, rt_tol_) in zip(lib, rt_tol)
        ischemicalequal(ion, getproperty(i, libchemical)) && continue
        irt = getrt(i)
        isnan(irt) && continue
        in(irt, rt_tol_) || continue
        imz1 = getmz(i)
        in(imz1, mz_tol) || continue
        push!(c, getproperty(i, libchemical))
    end
    c
end
isobars_rt_table(ion::AbstractChemical, lib::Table; fwhm_exp = 0.2, libid = true, libchemical = :AdductIon, libmz = :MZ, librt = :RT, libfwhm = :FWHM, mz_tol = crit(0.01, 20e-6)) = 
    isobars_rt_table(ion, rt(ion), mz(ion), lib; fwhm_exp, libid, libchemical = Symbol(libchemical), libmz = isnothing(libmz) ? nothing : Symbol(libmz), librt = isnothing(librt) ? nothing : Symbol(librt), libfwhm = isnothing(libfwhm) ? nothing : Symbol(libfwhm), mz_tol)
function isobars_rt_table(ion::AbstractChemical, rrt, rmz1, lib::Table; fwhm_exp = 0.2, libid = true, libchemical = :AdductIon, libmz = :MZ, librt = :RT, libfwhm = :FWHM, mz_tol = crit(0.01, 20e-6))
    isnan(rrt) && return nothing
    getmz = isnothing(libmz) ? (i -> mz(getproperty(i, libchemical))) : (i -> getproperty(i, libmz))
    getrt = isnothing(librt) ? (i -> rt(getproperty(i, libchemical))) : (i -> getproperty(i, librt))
    fwhm = isnothing(libfwhm) ? [fwhm_exp for _ in eachindex(lib)] : (getproperty(lib, libfwhm) .+ fwhm_exp) ./ 2
    mz_tol = union(makecrit_delta(crit(mz_tol), rmz1)...)
    rt_tol = [union(makecrit_delta(crit(x), rrt)...) for x in fwhm]
    lid = Int[]
    Δrt = Float64[]
    Δrtp = Float64[]
    Δmz = Float64[]
    Δmzp = Float64[]
    for (j, (i, rt_tol_)) in enumerate(zip(lib, rt_tol))
        ischemicalequal(ion, getproperty(i, libchemical)) && continue
        irt = getrt(i)
        isnan(irt) && continue
        in(irt, rt_tol_) || continue
        imz1 = getmz(i)
        in(imz1, mz_tol) || continue
        push!(lid, j)
        push!(Δrt, rrt - irt)
        push!(Δrtp, last(Δrt) / irt * 100)
        push!(Δmz, rmz1 - imz1)
        push!(Δmzp, last(Δmz) / imz1 * 1e6)
    end
    id = sortperm(Δrt; by = abs)
    lid = lid[id]
    Δrt = Δrt[id]
    Δrtp = Δrtp[id]
    Δmz = Δmz[id]
    Δmzp = Δmzp[id]
    libid ? Table(; Isobar = getproperty(lib, libchemical)[lid], ΔRT = Δrt, var"ΔRT(%)" = Δrtp, ΔMZ = Δmz, var"ΔMZ(ppm)" = Δmzp, LibID = lid) : 
        Table(; Isobar = getproperty(lib, libchemical)[lid], ΔRT = Δrt, var"ΔRT(%)" = Δrtp, ΔMZ = Δmz, var"ΔMZ(ppm)" = Δmzp)
end

function isobars_rt_table(exp::Vector, lib::Table; fwhm_exp = 0.2, libid = true, libchemical = :AdductIon, libmz = :MZ, librt = :RT, libfwhm = :FWHM, mz_tol = crit(0.01, 20e-6))
    fwhm_exp = collecfwhm(fwhm_exp, exp)
    libchemical = Symbol(libchemical)
    libmz = isnothing(libmz) ? nothing : Symbol(libmz)
    librt = isnothing(librt) ? nothing : Symbol(librt)
    libfwhm = isnothing(libfwhm) ? nothing : Symbol(libfwhm)
    if length(exp) < Threads.nthreads()
        mapreduce(vcat, exp, fwhm_exp) do x, y
            it = isobars_rt_table(x, lib; fwhm_exp = y, libid, libchemical, libmz, librt, libfwhm, mz_tol)
            Table((AdductIon = vectorize(x, length(it)), ), it)
        end |> Table
    else
        t = Vector{Table}(undef, length(exp))
        Threads.@threads for i in eachindex(t)
            it = isobars_rt_table(exp[i], lib; fwhm_exp = fwhm_exp[i], libid, libchemical, libmz, librt, libfwhm, mz_tol)
            t[i] = Table((AdductIon = vectorize(exp[i], length(it)), ), it)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end

function isobars_rt_table(exp::Table, lib::Table; expid = true, expchemical = :AdductIon, expmz = :MZ, exprt = :RT, expfwhm = :FWHM, libid = true, libchemical = :AdductIon, libmz = :MZ, librt = :RT, libfwhm = :FWHM, mz_tol = crit(0.01, 20e-6))
    libchemical = Symbol(libchemical)
    libmz = isnothing(libmz) ? nothing : Symbol(libmz)
    librt = isnothing(librt) ? nothing : Symbol(librt)
    libfwhm = isnothing(libfwhm) ? nothing : Symbol(libfwhm)
    expchemical = Symbol(expchemical)
    expmz = isnothing(expmz) ? nothing : Symbol(expmz)
    exprt = isnothing(exprt) ? nothing : Symbol(exprt)
    expfwhm = isnothing(expfwhm) ? nothing : Symbol(expfwhm)
    rrt = isnothing(exprt) ? rt.(getproperty(exp, expchemical)) : getproperty(exp, exprt)
    rmz1 = isnothing(expmz) ? mz.(getproperty(exp, expchemical)) : getproperty(exp, expmz)
    fwhm_exp = isnothing(expfwhm) ? [0.2 for _ in eachindex(exp)] : getproperty(exp, expfwhm)
    if length(exp) < Threads.nthreads()
        mapreduce(vcat, eachindex(exp), getproperty(exp, expchemical), fwhm_exp, rrt, rmz1) do i, x, y, xrt, xmz1
            it = isobars_rt_table(x, xrt, xmz1, lib; fwhm_exp = y, libid, libchemical, librt, libmz, libfwhm, mz_tol)
            Table(expid ? (ExpID = vectorize(i, length(it)), AdductIon = vectorize(x, length(it)), ) : (AdductIon = vectorize(x, length(it)), ), it)
        end |> Table
    else
        t = Vector{Table}(undef, length(exp))
        Threads.@threads for i in eachindex(t)
            it = isobars_rt_table(getproperty(exp, expchemical)[i], rrt[i], rmz1[i], lib; fwhm_exp = fwhm_exp[i], libid, libchemical, librt, libmz, libfwhm, mz_tol)
            t[i] = Table(expid ? (ExpID = vectorize(i, length(it)), AdductIon = vectorize(getproperty(exp, expchemical)[i], length(it)), ) : (AdductIon = vectorize(getproperty(exp, expchemical)[i], length(id)), ), it)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end