function getisobars(ion::AbstractIon, lib::AbstractVector{T}; rt_tol = 0.3, mz_tol = crit(0.01, 20e-6)) where {T <: Union{AbstractIon, IonCluster}}
    isnothing(rt(ion)) && return nothing
    rrt = rt(ion)
    rmz1 = mz(ion)
    mz_tol = union(crit(real_interval(mz_tol))(rmz1)...)
    rt_tol = union(crit(real_interval(rt_tol))(rrt)...)
    c = T[]
    for i in lib
        isionequal(ion, i) && continue
        irt = rt(i)
        isnothing(irt) && continue
        in(irt - rrt, rt_tol) || continue
        imz1 = mz(i)
        in(imz1 - rmz1, mz_tol) || continue
        push!(c, i)
    end
    c
end
isobar_rt_mz1(ion::AbstractIon, lib::AbstractVector{T}; rt_tol = 0.3, mz_tol = crit(0.01, 20e-6)) where {T <: Union{AbstractIon, IonCluster}} = 
    isobar_rt_mz1(ion, rt(ion), mz(ion), lib; rt_tol, mz_tol)
function isobar_rt_mz1(ion::AbstractIon, rrt, rmz1, lib::AbstractVector{T}; rt_tol = 0.3, mz_tol = crit(0.01, 20e-6)) where {T <: Union{AbstractIon, IonCluster}}
    isnothing(rrt) && return nothing
    mz_tol = union(crit(real_interval(mz_tol))(rmz1)...)
    rt_tol = union(crit(real_interval(rt_tol))(rrt)...)
    c = T[]
    Δrt = Float64[]
    Δmz = Float64[]
    for i in lib
        isionequal(ion, i) && continue
        irt = rt(i)
        isnothing(irt) && continue
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

function isobar_rt_mz1(exp::Vector, lib::AbstractVector{T}; kwargs...) where {T <: Union{AbstractIon, IonCluster}}
    if length(exp) < Threads.nthreads()
        mapreduce(vcat, exp) do x
            it = isobar_rt_mz1(x, lib; kwargs...)
            Table((Ion = vectorize(x, length(it)), ), it)
        end |> Table
    else
        t = Vector{Table}(undef, length(exp))
        Threads.@threads for i in eachindex(t)
            it = isobar_rt_mz1(exp[i], lib; kwargs...)
            t[i] = Table((Ion = vectorize(exp[i], length(it)), ), it)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end

function isobar_rt_mz1(exp::Table, lib::AbstractVector{T}; expid = true, expchemical = :Ion, expmz = :MZ, exprt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6)) where {T <: Union{AbstractIon, IonCluster}}
    expchemical = Symbol(expchemical)
    expmz = isnothing(expmz) ? nothing : Symbol(expmz)
    exprt = isnothing(exprt) ? nothing : Symbol(exprt)
    rrt = isnothing(exprt) ? rt.(getproperty(exp, expchemical)) : getproperty(exp, exprt)
    rmz1 = isnothing(expmz) ? mz.(getproperty(exp, expchemical)) : getproperty(exp, expmz)
    if length(exp) < Threads.nthreads()
        mapreduce(vcat, eachindex(exp), getproperty(exp, expchemical), rrt, rmz1) do i, x, xrt, xmz1
            it = isobar_rt_mz1(x, xrt, xmz1, lib; rt_tol, mz_tol)
            Table(expid ? (ExpID = vectorize(i, length(it)), Ion = vectorize(x, length(it)), ) : (Ion = vectorize(x, length(it)), ), it)
        end |> Table
    else
        t = Vector{Table}(undef, length(exp))
        Threads.@threads for i in eachindex(t)
            it = isobar_rt_mz1(getproperty(exp, expchemical)[i], rrt[i], rmz1[i], lib; rt_tol, mz_tol)
            t[i] = Table(expid ? (ExpID = vectorize(i, length(it)), Ion = vectorize(getproperty(exp, expchemical)[i], length(it)), ) : (Ion = vectorize(getproperty(exp, expchemical)[i], length(id)), ), it)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end

function getisobars(ion::AbstractIon, lib::Table; libchemical = :Ion, libmz = :MZ, librt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6))
    libchemical = Symbol(libchemical)
    libmz = isnothing(libmz) ? nothing : Symbol(libmz)
    librt = isnothing(librt) ? nothing : Symbol(librt)
    isnothing(rt(ion)) && return nothing
    getmz = isnothing(libmz) ? (i -> mz(getproperty(i, libchemical))) : (i -> getproperty(i, libmz))
    getrt = isnothing(librt) ? (i -> rt(getproperty(i, libchemical))) : (i -> getproperty(i, librt))
    rrt = rt(ion)
    rmz1 = mz(ion)
    mz_tol = union(crit(real_interval(mz_tol))(rmz1)...)
    rt_tol = union(crit(real_interval(rt_tol))(rrt)...)
    c = eltype(getproperty(lib, libchemical))[]
    for i in lib
        isionequal(ion, getproperty(i, libchemical)) && continue
        irt = getrt(i)
        isnothing(irt) && continue
        in(irt - rrt, rt_tol) || continue
        imz1 = getmz(i)
        in(imz1 - rmz1, mz_tol) || continue
        push!(c, i)
    end
    c
end
isobar_rt_mz1(ion::AbstractIon, lib::Table; libid = true, libchemical = :Ion, libmz = :MZ, librt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6)) = 
    isobar_rt_mz1(ion, rt(ion), mz(ion), lib; libid, libchemical = Symbol(libchemical), libmz = isnothing(libmz) ? nothing : Symbol(libmz), librt = isnothing(librt) ? nothing : Symbol(librt), rt_tol, mz_tol)
function isobar_rt_mz1(ion::AbstractIon, rrt, rmz1, lib::Table; libid = true, libchemical = :Ion, libmz = :MZ, librt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6))
    isnothing(rrt) && return nothing
    getmz = isnothing(libmz) ? (i -> mz(getproperty(i, libchemical))) : (i -> getproperty(i, libmz))
    getrt = isnothing(librt) ? (i -> rt(getproperty(i, libchemical))) : (i -> getproperty(i, librt))
    mz_tol = union(crit(real_interval(mz_tol))(rmz1)...)
    rt_tol = union(crit(real_interval(rt_tol))(rrt)...)
    lid = Int[]
    Δrt = Float64[]
    Δmz = Float64[]
    for (j, i) in enumerate(lib)
        isionequal(ion, getproperty(i, libchemical)) && continue
        irt = getrt(i)
        isnothing(irt) && continue
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

function isobar_rt_mz1(exp::Vector, lib::Table; libid = true, libchemical = :Ion, libmz = :MZ, librt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6))
    libchemical = Symbol(libchemical)
    libmz = isnothing(libmz) ? nothing : Symbol(libmz)
    librt = isnothing(librt) ? nothing : Symbol(librt)
    if length(v) < Threads.nthreads()
        mapreduce(vcat, v) do x, m
            it = isobar_rt_mz1(x, lib; lidid, libchemical, libmz, librt, rt_tol, mz_tol)
            Table((Ion = vectorize(x, length(it)), ), it)
        end |> Table
    else
        t = Vector{Table}(undef, length(v))
        Threads.@threads for i in eachindex(t)
            it = isobar_rt_mz1(v[i], lib; lidid, libchemical, libmz, librt, rt_tol, mz_tol)
            t[i] = Table((Ion = vectorize(v[i], length(it)), ), it)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end

function isobar_rt_mz1(exp::Table, lib::Table; expid = true, expchemical = :Ion, expmz = :MZ, exprt = :RT, libid = true, libchemical = :Ion, libmz = :MZ, librt = :RT, rt_tol = 0.3, mz_tol = crit(0.01, 20e-6))
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
            it = isobar_rt_mz1(x, xrt, xmz1, lib; libid, libchemical, librt, libmz, rt_tol, mz_tol)
            Table(expid ? (ExpID = vectorize(i, length(it)), Ion = vectorize(x, length(it)), ) : (Ion = vectorize(x, length(it)), ), it)
        end |> Table
    else
        t = Vector{Table}(undef, length(exp))
        Threads.@threads for i in eachindex(t)
            it = isobar_rt_mz1(getproperty(exp, expchemical)[i], rrt[i], rmz1[i], lib; libid, libchemical, librt, libmz, rt_tol, mz_tol)
            t[i] = Table(expid ? (ExpID = vectorize(i, length(it)), Ion = vectorize(getproperty(exp, expchemical)[i], length(it)), ) : (Ion = vectorize(getproperty(exp, expchemical)[i], length(id)), ), it)
        end
        Table(; (map(propertynames(t[1])) do p
            p => ChainedVector(getproperty.(t, p))
        end)...)
    end
end