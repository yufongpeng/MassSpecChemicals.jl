"""
    getisobars_rt(ion::AbstractIon, lib::AbstractVector{T}; rt_tol = 0.3, mz_tol = crit(0.01, 20e-6))
    getisobars_rt(ion::AbstractIon, lib::Table; libchemical = :Ion, libmz = :MZ, librt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6))

Return a subvector of `lib` such that all elements are isobaric to `ion` within specific tolerance of rt and mz. 

Keyword arguments `libchemical`, `libmz`, and `librt` specify the columns of `Ion` objects, m/z values, and rt values. 

Tolerance can be a number or criteria, representing maximum allowed difference. For criteria with multiple allowed intervals, lying in any of them is considered to fulfill the criteria. 
"""
function getisobars_rt(ion::AbstractIon, lib::AbstractVector{T}; rt_tol = 0.3, mz_tol = crit(0.01, 20e-6)) where {S <: AbstractIon, T <: Union{AbstractIon, Isobars{S}}}
    isnan(rt(ion)) && return T[]
    rrt = rt(ion)
    rmz1 = mz(ion)
    mz_tol = union(crit(real_interval(mz_tol))(rmz1)...)
    rt_tol = union(crit(real_interval(rt_tol))(rrt)...)
    c = T[]
    for i in lib
        ischemicalequal(ion, i) && continue
        irt = rt(i)
        isnan(irt) && continue
        in(irt - rrt, rt_tol) || continue
        imz1 = mz(i)
        in(imz1 - rmz1, mz_tol) || continue
        push!(c, i)
    end
    c
end

"""
    isobartable_rt(ion::AbstractIon, lib::AbstractVector{T}; rt_tol = 0.3, mz_tol = crit(0.01, 20e-6))
    isobartable_rt(exp::Vector, lib::AbstractVector{T}; kwargs...)
    isobartable_rt(exp::Table, lib::AbstractVector{T}; expid = true, expchemical = :Ion, expmz = :MZ, exprt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6))
    isobartable_rt(ion::AbstractIon, lib::Table; libid = true, libchemical = :Ion, libmz = :MZ, librt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6))
    isobartable_rt(exp::Vector, lib::Table; libid = true, libchemical = :Ion, libmz = :MZ, librt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6))
    isobartable_rt(exp::Table, lib::Table; expid = true, expchemical = :Ion, expmz = :MZ, exprt = :RT, libid = true, libchemical = :Ion, libmz = :MZ, librt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6))

Return a `Table` such that all rows representing an ion from `lib` which are isobaric to `ion` or ions from `exp` within specific tolerance of rt and mz. 

Tolerance can be a number or criteria, representing maximum allowed difference. For criteria with multiple allowed intervals, lying in any of them is considered to fulfill the criteria. 

Keyword arguments `libchemical`, `libmz`, and `librt` specify the columns of `Ion` objects, m/z values, and rt values from `lib`. 

Keyword arguments `expchemical`, `expmz`, and `exprt` specify the columns of `Ion` objects, m/z values, and rt values from `exp`. 

Keyword arguments `libid` and `expid` determine whether creates columns to store row id of `lib` or `exp`.
"""
isobartable_rt(ion::AbstractIon, lib::AbstractVector{T}; rt_tol = 0.3, mz_tol = crit(0.01, 20e-6)) where {S <: AbstractIon, T <: Union{AbstractIon, Isobars{S}}} = 
    isobartable_rt(ion, rt(ion), mz(ion), lib; rt_tol, mz_tol)
function isobartable_rt(ion::AbstractIon, rrt, rmz1, lib::AbstractVector{T}; rt_tol = 0.3, mz_tol = crit(0.01, 20e-6)) where {S <: AbstractIon, T <: Union{AbstractIon, Isobars{S}}}
    isnan(rrt) && return nothing
    mz_tol = union(crit(real_interval(mz_tol))(rmz1)...)
    rt_tol = union(crit(real_interval(rt_tol))(rrt)...)
    c = T[]
    Δrt = Float64[]
    Δmz = Float64[]
    for i in lib
        ischemicalequal(ion, i) && continue
        irt = rt(i)
        isnan(irt) && continue
        drt = irt - rrt
        in(drt, rt_tol) || continue
        imz1 = mz(i)
        dmz = imz1 - rmz1
        in(dmz, mz_tol) || continue
        push!(c, i)
        push!(Δrt, drt)
        push!(Δmz, dmz)
    end
    id = sortperm(Δrt)
    c = c[id]
    Δrt = Δrt[id]
    Δmz = Δmz[id]
    Table(; Isobar = c, ΔRT = Δrt, var"ΔRT(%)" = Δrt ./ rrt .* 100, ΔMZ = Δmz, var"ΔMZ(ppm)" = Δmz ./ rmz1 .* 1e6)
end

function isobartable_rt(exp::Vector, lib::AbstractVector{T}; kwargs...) where {S <: AbstractIon, T <: Union{AbstractIon, Isobars{S}}}
    if length(exp) < Threads.nthreads()
        mapreduce(vcat, exp) do x
            it = isobartable_rt(x, lib; kwargs...)
            Table((Ion = vectorize(x, length(it)), ), it)
        end |> Table
    else
        t = Vector{Table}(undef, length(exp))
        Threads.@threads for i in eachindex(t)
            it = isobartable_rt(exp[i], lib; kwargs...)
            t[i] = Table((Ion = vectorize(exp[i], length(it)), ), it)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end

function isobartable_rt(exp::Table, lib::AbstractVector{T}; expid = true, expchemical = :Ion, expmz = :MZ, exprt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6)) where {S <: AbstractIon, T <: Union{AbstractIon, Isobars{S}}}
    expchemical = Symbol(expchemical)
    expmz = isnothing(expmz) ? nothing : Symbol(expmz)
    exprt = isnothing(exprt) ? nothing : Symbol(exprt)
    rrt = isnothing(exprt) ? rt.(getproperty(exp, expchemical)) : getproperty(exp, exprt)
    rmz1 = isnothing(expmz) ? mz.(getproperty(exp, expchemical)) : getproperty(exp, expmz)
    if length(exp) < Threads.nthreads()
        mapreduce(vcat, eachindex(exp), getproperty(exp, expchemical), rrt, rmz1) do i, x, xrt, xmz1
            it = isobartable_rt(x, xrt, xmz1, lib; rt_tol, mz_tol)
            Table(expid ? (ExpID = vectorize(i, length(it)), Ion = vectorize(x, length(it)), ) : (Ion = vectorize(x, length(it)), ), it)
        end |> Table
    else
        t = Vector{Table}(undef, length(exp))
        Threads.@threads for i in eachindex(t)
            it = isobartable_rt(getproperty(exp, expchemical)[i], rrt[i], rmz1[i], lib; rt_tol, mz_tol)
            t[i] = Table(expid ? (ExpID = vectorize(i, length(it)), Ion = vectorize(getproperty(exp, expchemical)[i], length(it)), ) : (Ion = vectorize(getproperty(exp, expchemical)[i], length(id)), ), it)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end

function getisobars_rt(ion::AbstractIon, lib::Table; libchemical = :Ion, libmz = :MZ, librt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6))
    libchemical = Symbol(libchemical)
    libmz = isnothing(libmz) ? nothing : Symbol(libmz)
    librt = isnothing(librt) ? nothing : Symbol(librt)
    isnan(rt(ion)) && return eltype(getproperty(lib, libchemical))[]
    getmz = isnothing(libmz) ? (i -> mz(getproperty(i, libchemical))) : (i -> getproperty(i, libmz))
    getrt = isnothing(librt) ? (i -> rt(getproperty(i, libchemical))) : (i -> getproperty(i, librt))
    rrt = rt(ion)
    rmz1 = mz(ion)
    mz_tol = union(crit(real_interval(mz_tol))(rmz1)...)
    rt_tol = union(crit(real_interval(rt_tol))(rrt)...)
    c = eltype(getproperty(lib, libchemical))[]
    for i in lib
        ischemicalequal(ion, getproperty(i, libchemical)) && continue
        irt = getrt(i)
        isnan(irt) && continue
        in(irt - rrt, rt_tol) || continue
        imz1 = getmz(i)
        in(imz1 - rmz1, mz_tol) || continue
        push!(c, i)
    end
    c
end
isobartable_rt(ion::AbstractIon, lib::Table; libid = true, libchemical = :Ion, libmz = :MZ, librt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6)) = 
    isobartable_rt(ion, rt(ion), mz(ion), lib; libid, libchemical = Symbol(libchemical), libmz = isnothing(libmz) ? nothing : Symbol(libmz), librt = isnothing(librt) ? nothing : Symbol(librt), rt_tol, mz_tol)
function isobartable_rt(ion::AbstractIon, rrt, rmz1, lib::Table; libid = true, libchemical = :Ion, libmz = :MZ, librt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6))
    isnan(rrt) && return nothing
    getmz = isnothing(libmz) ? (i -> mz(getproperty(i, libchemical))) : (i -> getproperty(i, libmz))
    getrt = isnothing(librt) ? (i -> rt(getproperty(i, libchemical))) : (i -> getproperty(i, librt))
    mz_tol = union(crit(real_interval(mz_tol))(rmz1)...)
    rt_tol = union(crit(real_interval(rt_tol))(rrt)...)
    lid = Int[]
    Δrt = Float64[]
    Δmz = Float64[]
    for (j, i) in enumerate(lib)
        ischemicalequal(ion, getproperty(i, libchemical)) && continue
        irt = getrt(i)
        isnan(irt) && continue
        drt = irt - rrt
        in(drt, rt_tol) || continue
        imz1 = getmz(i)
        dmz = imz1 - rmz1
        in(dmz, mz_tol) || continue
        push!(lid, j)
        push!(Δrt, drt)
        push!(Δmz, dmz)
    end
    id = sortperm(Δrt)
    lid = lid[id]
    Δrt = Δrt[id]
    Δmz = Δmz[id]
    libid ? Table(; Isobar = getproperty(lib, libchemical)[lid], ΔRT = Δrt, var"ΔRT(%)" = Δrt ./ rrt .* 100, ΔMZ = Δmz, var"ΔMZ(ppm)" = Δmz ./ rmz1 .* 1e6, LibID = lid) : 
        Table(; Isobar = getproperty(lib, libchemical)[lid], ΔRT = Δrt, var"ΔRT(%)" = Δrt ./ rrt .* 100, ΔMZ = Δmz, var"ΔMZ(ppm)" = Δmz ./ rmz1 .* 1e6)
end

function isobartable_rt(exp::Vector, lib::Table; libid = true, libchemical = :Ion, libmz = :MZ, librt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6))
    libchemical = Symbol(libchemical)
    libmz = isnothing(libmz) ? nothing : Symbol(libmz)
    librt = isnothing(librt) ? nothing : Symbol(librt)
    if length(v) < Threads.nthreads()
        mapreduce(vcat, exp) do x, m
            it = isobartable_rt(x, lib; libid, libchemical, libmz, librt, rt_tol, mz_tol)
            Table((Ion = vectorize(x, length(it)), ), it)
        end |> Table
    else
        t = Vector{Table}(undef, length(exp))
        Threads.@threads for i in eachindex(t)
            it = isobartable_rt(exp[i], lib; libid, libchemical, libmz, librt, rt_tol, mz_tol)
            t[i] = Table((Ion = vectorize(exp[i], length(it)), ), it)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end

function isobartable_rt(exp::Table, lib::Table; expid = true, expchemical = :Ion, expmz = :MZ, exprt = :RT, libid = true, libchemical = :Ion, libmz = :MZ, librt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6))
    libchemical = Symbol(libchemical)
    libmz = isnothing(libmz) ? nothing : Symbol(libmz)
    librt = isnothing(librt) ? nothing : Symbol(librt)
    expchemical = Symbol(expchemical)
    expmz = isnothing(expmz) ? nothing : Symbol(expmz)
    exprt = isnothing(exprt) ? nothing : Symbol(exprt)
    rrt = isnothing(exprt) ? rt.(getproperty(exp, expchemical)) : getproperty(exp, exprt)
    rmz1 = isnothing(expmz) ? mz.(getproperty(exp, expchemical)) : getproperty(exp, expmz)
    if length(exp) < Threads.nthreads()
        mapreduce(vcat, eachindex(exp), getproperty(exp, expchemical), rrt, rmz1) do i, x, xrt, xmz1
            it = isobartable_rt(x, xrt, xmz1, lib; libid, libchemical, librt, libmz, rt_tol, mz_tol)
            Table(expid ? (ExpID = vectorize(i, length(it)), Ion = vectorize(x, length(it)), ) : (Ion = vectorize(x, length(it)), ), it)
        end |> Table
    else
        t = Vector{Table}(undef, length(exp))
        Threads.@threads for i in eachindex(t)
            it = isobartable_rt(getproperty(exp, expchemical)[i], rrt[i], rmz1[i], lib; libid, libchemical, librt, libmz, rt_tol, mz_tol)
            t[i] = Table(expid ? (ExpID = vectorize(i, length(it)), Ion = vectorize(getproperty(exp, expchemical)[i], length(it)), ) : (Ion = vectorize(getproperty(exp, expchemical)[i], length(id)), ), it)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end