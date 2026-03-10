"""
    fwhm_mz(ms::AbstractMSAnalyzer, mz::Real)

FWHM of `mz`.
"""
fwhm_mz(ms::MSAnalyzer, mz) = ms.fwhm_mz(mz)
fwhm_mz(ms::Quadrupole, mz) = ms.fwhm 
fwhm_mz(ms::QuadrupoleIonTrap, mz) = ms.fwhm 
fwhm_mz(ms::LinearIonTrap, mz) = ms.fwhm 
fwhm_mz(ms::TOF, mz) = sqrt(mz * (3 * ms.mz50 + mz)) / ms.resolution
fwhm_mz(ms::Orbitrap, mz) = mz / ms.resolution * sqrt((3 * mz + ms.mz50 - 800) / (ms.mz50 - 200))
fwhm_mz(ms::FTICR, mz) = mz / ms.resolution * (mz + ms.mz50 - 400) / (ms.mz50 - 200)

"""
    resolving_power(ms::AbstractMSAnalyzer, mz::Real)

Resolving power of `ms` at `mz`.
"""
resolving_power(ms::AbstractMSAnalyzer, mz) = mz / fwhm_mz(ms, mz)

gaussian_norm(x::Real, μ::Real, σ2::Real) = exp(-(x - μ) ^ 2 / 2 / σ2)
supergaussian_norm(x::Real, μ::Real, σ::Real, n::Real) = exp(-((x - μ) / sqrt(2) / σ) ^ n)
sgaussian_norm(ν::Real, t) = exp(-ν ^ 2 / 2 / t)
sgaussian_norm(ν::AbstractRange, t::T) where T = sgaussian_norm!(zeros(T, length(ν)), ν, t)
function sgaussian_norm!(u::AbstractVector, ν::AbstractRange, t) 
    delta = step(ν)
    start = ν.start
    num = exp(-start ^ 2 / 2t)
    deltanum1 = exp(-delta ^ 2 / 2t)
    deltanum2 = exp(-delta / t)
    prev_ν = start - delta
    prev_num = num / deltanum1 / deltanum2 ^ prev_ν
    for i in eachindex(ν)
        u[i] = prev_num * deltanum1 * deltanum2 ^ prev_ν
        prev_num = u[i]
        prev_ν += delta
    end
    u
end

beseel_max(::Type{Float64}) = 713.0
beseel_max(::Type{Float32}) = Float32(91)
dgaussian_norm(ν::Real, t::T) where T = (t = min(t, beseel_max(T)); exp(-t) * besseli(ν, t))
dgaussian_norm(ν::AbstractRange, t::T) where T = dgaussian!(zeros(T, length(ν)), ν, t)
function dgaussian_norm!(u::AbstractVector, ν::AbstractRange, t) 
    t = min(t, beseel_max(T))
    u = besseli!(u, ν, t)
    f = exp(-t)
    u .*= f 
end

(::DiscreteGaussianWindow)(ν, t) = dgaussian_norm(ν, t)
(::SampledGaussianWindow)(ν, t) = sgaussian_norm(ν, t)
(::GaussianWindow)(x, μ, σ2) = gaussian_norm(x, μ, σ2)
(f::SuperGaussianWindow)(x, μ, σ) = supergaussian_norm(x, μ, σ, f.power)
function (f::FixedTaperTukeyWindow)(x, μ, fwhm)
    taperwidth = min(f.taperwidth, fwhm)
    p = abs(x - μ)
    plateau = fwhm - taperwidth
    if p > plateau / 2 + taperwidth
        0
    elseif p > plateau / 2
        (cos(p * π / taperwidth) + 1) / 2
    else
        1
    end
end

function (f::TukeyWindow)(x, μ, fwhm)
    width = fwhm / (2 - f.taperproportion)
    plateau = width * (1 - f.taperproportion)
    p = abs(x - μ) 
    if p > width
        0
    elseif p > plateau
        taperwidth = width - plateau
        p = p - plateau
        (cos(p * π / taperwidth) + 1) / 2
    else
        1
    end
end

(::CosineWindow)(x, μ, fwhm) = abs(x - μ) > fwhm ? 0 : (cos(π * abs(x - μ) / fwhm) + 1) / 2
(::RectWindow)(x, μ, fwhm) = abs(x - μ) < fwhm / 2 ? 1 : 0
function (f::PowerCosineWindow)(x, μ, fwhm)
    p = (x - μ) * 2 / fwhm * acos(2 ^ (1 - 1 / f.power) - 1)
    # p = (x - μ) * 2 / fwhm
    abs(p) > π ? 0 : (0.5 + cos(p) / 2) ^ f.power
end
    
function (f::GaussianTailedUniformWindow)(x, μ, fwhm)
    p = abs(x - μ) 
    p == 0 && return 1
    taperwidth = min(f.taperwidth, fwhm)
    plateau = fwhm - taperwidth
    if plateau == 0
        gaussian_norm(p, 0, fwhm ^ 2 / log(256))
    elseif p > plateau / 2
        gaussian_norm(p - plateau / 2, 0, taperwidth ^ 2 / log(256))
    else
        1
    end
end

function (f::GaussianDilatedWindow)(x, μ, fwhm)
    p = abs(x - μ) 
    p == 0 && return 1
    fwhm <= f.min_fwhm && return sgaussian_norm(p, fwhm ^ 2 / log(256))
    a = (fwhm - f.min_fwhm) / 2
    w = (fwhm - f.taperwidth) / 2
    t = p / w
    if t > 1
        sgaussian_norm(p - a, f.min_fwhm ^ 2 / log(256))
    elseif a > 0.75 * w && p ^ 2 <= (3w - 4a) / (4a - 2w)
        1
    else
        sgaussian_norm((a - w / 2) * t ^ 4 + (1.5 * w - 2 * a) * t ^ 2, f.min_fwhm ^ 2 / log(256))
    end
end

window_parameter(msanalyzer::AbstractMSAnalyzer{W}, μ) where {W <: GaussianTailedUniformWindow} = (μ, fwhm_mz(msanalyzer, μ))
window_parameter(msanalyzer::AbstractMSAnalyzer{W}, μ) where {W <: GaussianWindow} = (μ, fwhm_mz(msanalyzer, μ) ^ 2 / log(256))
window_parameter(msanalyzer::AbstractMSAnalyzer{W}, μ) where {W <: FixedTaperTukeyWindow} = (μ, fwhm_mz(msanalyzer, μ))
window_parameter(msanalyzer::AbstractMSAnalyzer{W}, μ) where {W <: TukeyWindow} = (μ, fwhm_mz(msanalyzer, μ))
window_parameter(msanalyzer::AbstractMSAnalyzer{W}, μ) where {W <: CosineWindow} = (μ, fwhm_mz(msanalyzer, μ))
window_parameter(msanalyzer::AbstractMSAnalyzer{W}, μ) where {W <: PowerCosineWindow} = (μ, fwhm_mz(msanalyzer, μ))
window_parameter(msanalyzer::AbstractMSAnalyzer{W}, μ) where {W <: RectWindow} = (μ, fwhm_mz(msanalyzer, μ))
window_parameter(msanalyzer::AbstractMSAnalyzer{W}, μ) where {W <: SuperGaussianWindow} = (μ, fwhm_mz(msanalyzer, μ) / (sqrt(8) * log(2) ^ (1 / msanalyzer.window.power)))

function discrete_window(window::AbstractWindow, μ::Vector{T}, msanalyzer, binsize, nbin_multiplier, height) where T
    fwhm_ref = fwhm_mz(msanalyzer, first(μ))
    window_ref = discrete_window(window, fwhm_ref, binsize, nbin_multiplier, height)
    windows = [Vector{T}() for _ in μ]
    for (i, m) in enumerate(μ)
        fwhm = fwhm_mz(msanalyzer, m)
        if fwhm - fwhm_ref < binsize * nbin_multiplier
            windows[i] = window_ref 
        else
            windows[i] = discrete_window(window, fwhm, binsize, nbin_multiplier, height)
            window_ref = windows[i]
            fwhm_ref = fwhm
        end
    end
    windows
end

function discrete_window(window::CosineWindow, fwhm, binsize, nbin_multiplier, height)
    dfwtm = ceil(fwhm * 2 / binsize)
    p, r = divrem(dfwtm, nbin_multiplier)
    dfwtm = r > 0 ? (p + 1) * nbin_multiplier : dfwtm
    k = [window(x * binsize, 0, fwhm) for x in 0:dfwtm]
    vcat(reverse(k), k[begin + 1:end])
end

function discrete_window(window::PowerCosineWindow, fwhm, binsize, nbin_multiplier, height)
    dfwtm = ceil(fwhm * 2 / binsize)
    p, r = divrem(dfwtm, nbin_multiplier)
    dfwtm = r > 0 ? (p + 1) * nbin_multiplier : dfwtm
    k = [window(x * binsize, 0, fwhm) for x in 0:dfwtm]
    vcat(reverse(k), k[begin + 1:end])
end

function discrete_window(window::TukeyWindow, fwhm, binsize, nbin_multiplier, height)
    dfwtm = ceil(fwhm / (2 - window.taperproportion) / binsize)
    plateau = ceil(dfwtm * (1 - window.taperproportion))
    cosfwhm = dfwtm - plateau
    p, r = divrem(dfwtm, nbin_multiplier)
    dfwtm = r > 0 ? (p + 1) * nbin_multiplier : dfwtm
    k = [CosineWindow()(x, 0, cosfwhm) for x in 0:dfwtm - plateau]
    vcat(reverse(k), [1 for _ in 1:(2 * plateau)], k[begin + 1:end])
end

function discrete_window(window::RectWindow, fwhm, binsize, nbin_multiplier, height)
    width = ceil(fwhm / 2 / binsize)
    p, r = divrem(width, nbin_multiplier)
    dfwtm = r > 0 ? (p + 1) * nbin_multiplier : width
    k = [Int(x <= width) for x in 0:dfwtm]
    vcat(reverse(k), k[begin + 1:end])
end

function discrete_window(window::FixedTaperTukeyWindow, fwhm, binsize, nbin_multiplier, height)
    dfwtm = ceil((fwhm + window.taperwidth) / 2 / binsize)
    p, r = divrem(dfwtm, nbin_multiplier)
    dfwtm = r > 0 ? (p + 1) * nbin_multiplier : dfwtm
    plateau = fwhm - window.taperwidth
    k = [window(x * binsize, 0, plateau) for x in 0:dfwtm]
    vcat(reverse(k), k[begin + 1:end])
end

function discrete_window(window::GaussianTailedUniformWindow, fwhm, binsize, nbin_multiplier, height)
    taperstart = ceil(Int, (fwhm - window.taperwidth) / 2 / binsize)
    offset = rem(taperstart, nbin_multiplier)
    fwhm = fwhm - taperstart * binsize * 2
    k = _discrete_window(SampledGaussianWindow(), fwhm, binsize, nbin_multiplier, height; offset)
    vcat(reverse(k), ones(2 * taperstart), k[begin + 1:end])
end

function discrete_window(window::GaussianDilatedWindow, fwhm, binsize, nbin_multiplier, height)
    fwhm <= window.min_fwhm && return discrete_window(SampledGaussianWindow(), fwhm, binsize, nbin_multiplier, height)
    taperstart = ceil(Int, (fwhm - window.taperwidth) / 2 / binsize)
    offset = rem(taperstart, nbin_multiplier)
    hwhm = window.min_fwhm / 2
    bsqrtt = hwhm / binsize
    t = bsqrtt ^ 2 / log(4)
    ihwhml = floor(Int, hwhm / binsize)
    ihwhmr = ceil(Int, hwhm / binsize)
    ihwhmr = ihwhmr > ihwhml ? ihwhmr : ihwhml + 1
    t = estimate_t(SampledGaussianWindow(), bsqrtt, ihwhml, ihwhmr)
    dfwtm = round(Int, sqrt(2t * log(1 / height))) + 1
    p, r = divrem(dfwtm, nbin_multiplier)
    dfwtm = r == offset ? dfwtm : r > offset ? (p + 1) * nbin_multiplier + offset : p * nbin_multiplier + offset
    if taperstart > 1
        new_min_fwhm = sqrt(t * log(256)) * binsize # fwhm - 2a
        a = ceil(Int, (fwhm - new_min_fwhm) / 2 / binsize)
        p = (0:(taperstart - 1)) ./ taperstart
        f1 = a - taperstart / 2
        f2 = 1.5 * taperstart - 2 * a
        f = @. sgaussian_norm(f1 * p ^ 4 + f2 * p ^ 2, t)
    else
        a = 0 
        f = [1]
    end
    k = SampledGaussianWindow()((taperstart - a):dfwtm + taperstart - a, t)
    k = reverse(vcat(f, k))
    m = 0
    for i in eachindex(k)
        if k[i] >= m 
            m = k[i]
        else
            k[i] = m
        end
    end
    vcat(k, reverse(k[begin:end - 1]))
end

discrete_window(window::GaussianWindow, fwhm, binsize, nbin_multiplier, height) = discrete_window(SampledGaussianWindow(), fwhm, binsize, nbin_multiplier, height)
function discrete_window(window::SuperGaussianWindow, fwhm, binsize, nbin_multiplier, height)
    hwhm = fwhm / 2
    bsqrtt = hwhm / binsize
    ihwhml = floor(Int, hwhm / binsize)
    ihwhmr = ceil(Int, hwhm / binsize)
    ihwhmr = ihwhmr > ihwhml ? ihwhmr : ihwhml + 1
    t = estimate_t(window, bsqrtt, ihwhml, ihwhmr)
    dfwtm = round(Int, sqrt(2) * t * log(1 / height) ^ (1 / window.power)) + 1
    p, r = divrem(dfwtm, nbin_multiplier)
    dfwtm = r > 0 ? (p + 1) * nbin_multiplier : dfwtm
    k = [window(x, 0, t) for x in 0:dfwtm]
    vcat(reverse(k), k[begin + 1:end])
end
function _discrete_window(window::SampledGaussianWindow, fwhm, binsize, nbin_multiplier, height; offset = 0)
    hwhm = fwhm / 2
    bsqrtt = hwhm / binsize
    ihwhml = floor(Int, hwhm / binsize)
    ihwhmr = ceil(Int, hwhm / binsize)
    ihwhmr = ihwhmr > ihwhml ? ihwhmr : ihwhml + 1
    t = estimate_t(window, bsqrtt, ihwhml, ihwhmr)
    dfwtm = round(Int, sqrt(2t * log(1 / height))) + 1
    p, r = divrem(dfwtm, nbin_multiplier)
    dfwtm = r == offset ? dfwtm : r > offset ? (p + 1) * nbin_multiplier + offset : p * nbin_multiplier + offset
    window(0:dfwtm, t)
end

function discrete_window(window::SampledGaussianWindow, fwhm, binsize, nbin_multiplier, height)
    k = _discrete_window(window, fwhm, binsize, nbin_multiplier, height)
    vcat(reverse(k), k[begin + 1:end])
end

function discrete_window(window::DiscreteGaussianWindow, fwhm::T, binsize, nbin_multiplier, height) where T 
    hwhm = fwhm / 2
    bsqrtt = hwhm / binsize
    t = (bsqrtt / log(4)) ^ 2
    if t > bessel_max(T)
        delta = sqrt(t / beseel_max(T))
        dfwtm = sqrt(2t * log(1 / height)) + delta
        p, r = divrem(dfwtm, nbin_multiplier * delta)
        dfwtm = r > 0 ? (p + 1) * nbin_multiplier * delta - r : dfwtm
        t = beseel_max(T)
        k = map(x -> window(x, t), 0:delta:dfwtm)
    else
        ihwhml = floor(Int, hwhm / binsize)
        ihwhmr = ceil(Int, hwhm / binsize)
        ihwhmr = ihwhmr > ihwhml ? ihwhmr : ihwhml + 1
        t = estimate_t(window, bsqrtt, ihwhml, ihwhmr)
        dfwtm = round(Int, sqrt(2t * log(1 / height))) + 1
        p, r = divrem(dfwtm, nbin_multiplier)
        dfwtm = r > 0 ? (p + 1) * nbin_multiplier : dfwtm
        k = window(0:dfwtm, t)
    end
    vcat(reverse(k), k[begin + 1:end])
end

function estimate_t(window::SampledGaussianWindow, bsqrtt, ihwhml, ihwhmr, log_n = 0)
    t = bsqrtt ^ 2 / log(4)
    log_n > 1000 && return t
    # p = window(0, t) / 2
    hl = window(ihwhml, t)
    hr = window(ihwhmr, t)
    if hr <= 0.5 <= hl 
        return t
    elseif 0.5 > hl
        bsqrtt += (0.5 - (hl + hr) / 2) / (hl - hr)
        return estimate_t(window, bsqrtt, ihwhml, ihwhmr, log_n + 1)
    else
        bsqrtt += (0.5 - (hl + hr) / 2) / (hl - hr)
        return estimate_t(window, bsqrtt, ihwhml, ihwhmr, log_n + 1)
    end
end

function estimate_t(window::SuperGaussianWindow, bsqrtt, ihwhml, ihwhmr, log_n = 0)
    t = bsqrtt / sqrt(2) / log(2) ^ (1 / window.power)
    log_n > 1000 && return t
    # p = window(0, 0, t) / 2
    hl = window(ihwhml, 0, t)
    hr = window(ihwhmr, 0, t)
    if hr <= 0.5 <= hl 
        return t
    elseif 0.5 > hl
        bsqrtt += (0.5 - (hl + hr) / 2) / (hl - hr)
        return estimate_t(window, bsqrtt, ihwhml, ihwhmr, log_n + 1)
    else
        bsqrtt += (0.5 - (hl + hr) / 2) / (hl - hr)
        return estimate_t(window, bsqrtt, ihwhml, ihwhmr, log_n + 1)
    end
end