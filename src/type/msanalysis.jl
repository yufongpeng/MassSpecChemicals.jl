"""
    Spectrum{T}

Discrete mass spectrum.

# Fields
* `spectrum::Vector{T}`: intensity values on the discretized m/z grid.
* `initial_mass::T`: the m/z value at the start of the grid.
* `binsize::T`: the m/z difference between adjacent grid points.
* `stepsize::Int`: number of grid bins represented by each point in the original spectrum.
* `table::Table`: source table of the ions.
"""
struct Spectrum{T}
    spectrum::Vector{T}
    initial_mass::T 
    binsize::T
    stepsize::Int
    table::Table
end

# Record procedure
# struct MSAnalysis 
#     procedure::Vector
#     table::Table
#     tables::Vector
# end

abstract type AbstractWindow end
abstract type ZeroCenteredWindow end

"""
    GaussianWindow <: AbstractWindow
    (::GaussianWindow)(x, μ, σ2)

Gaussian window evaluated at `x` with mean (center) `μ` and variance `σ2`.
"""
struct GaussianWindow <: AbstractWindow end

"""
    GaussianTailedUniformWindow{T} <: AbstractWindow
    (::GaussianTailedUniformWindow)(x, μ, fwhm)

Uniform window with Gaussian tails, evaluated at `x` with center `μ` and full width at half maximum `fwhm`.

# Fields
* `taperwidth`: width of each Gaussian tail.
"""
struct GaussianTailedUniformWindow{T} <: AbstractWindow
    taperwidth::T
end

"""
    FixedTaperTukeyWindow{T} <: AbstractWindow
    (::FixedTaperTukeyWindow)(x, μ, fwhm)

Tukey window with a fixed taper width, evaluated at `x` with center `μ` and full width at half maximum `fwhm`.

# Fields
* `taperwidth`: width of the tapered region.
"""
struct FixedTaperTukeyWindow{T} <: AbstractWindow 
    taperwidth::T
end

"""
    TukeyWindow{T} <: AbstractWindow
    (::TukeyWindow)(x, μ, fwhm)

Tukey window with a proportional taper, evaluated at `x` with center `μ` and full width at half maximum `fwhm`.

# Fields
* `taperproportion`: proportion of the window occupied by the tapered region.
"""
struct TukeyWindow{T} <: AbstractWindow 
    taperproportion::T
end

"""
    CosineWindow <: AbstractWindow
    (::CosineWindow)(x, μ, fwhm)

Cosine window evaluated at `x` with mean (center) `μ` and full width at half maximum `fwhm`.
"""
struct CosineWindow <: AbstractWindow end

"""
    RectWindow <: AbstractWindow
    (::RectWindow)(x, μ, fwhm)

Rectangular window evaluated at `x` with mean (center) `μ` and full width at half maximum `fwhm`.
"""
struct RectWindow <: AbstractWindow end

"""
    SampledGaussianWindow <: ZeroCenteredWindow
    (::SampledGaussianWindow)(x, t)

Gaussian window evaluated at `x` on a grid centered at zero, with variance `t`.
"""
struct SampledGaussianWindow <: ZeroCenteredWindow end

"""
    SuperGaussianWindow{T} <: AbstractWindow
    (::SuperGaussianWindow)(x, μ, σ)

Super-Gaussian window evaluated at `x` with center `μ` and scale `σ`. A power of 2 gives a standard Gaussian shape.

# Fields
* `power`: exponent that controls the shape of the window.
"""
struct SuperGaussianWindow{T} <: AbstractWindow 
    power::T
end

"""
    PowerCosineWindow{T} <: AbstractWindow
    (::PowerCosineWindow)(x, μ, fwhm)

Power-cosine (top-half cosine) window evaluated at `x` with center `μ` and full width at half maximum `fwhm`.

# Fields
* `power`: exponent applied to the cosine function.
"""
struct PowerCosineWindow{T} <: AbstractWindow 
    power::T
end

abstract type AbstractMSAnalyzer{W, M} end 

"""
    MSAnalyzer{W, M} <: AbstractMSAnalyzer{W, M} 

Generic mass-spectrometer analyzer.

# Fields
* `window::AbstractWindow`: window function used to model the analyzer response.
* `mz`: analyzed m/z value or values:
    * `Real`: one center m/z value.
    * `Vector{<:Real}`: multiple center m/z values.
    * `Tuple{<:Real, <:Real}`: an m/z range.
    * `Nothing`: no m/z value is specified.
* `accuracy`: m/z accuracy, specified as a number or a [`Criteria`](@ref).
* `fwhm_mz`: function that maps an m/z value to its full width at half maximum.
* `digits`: number of decimal places retained for m/z values. `nothing` means there is no limit.
* `stepsize`: step size for a discrete mass scan, such as a quadrupole scan. It must be a multiple of `10 ^ (-digits)`. If `digits` is `nothing`, `stepsize` must also be `nothing`.
"""
struct MSAnalyzer{W, M} <: AbstractMSAnalyzer{W, M}
    window::W
    mz::M
    accuracy
    fwhm_mz 
    digits
    stepsize
end

"""
    Quadrupole{W, M} <: AbstractMSAnalyzer{W, M} 

Quadrupole mass analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value or values:
    * `Real`: one center m/z value.
    * `Vector{<:Real}`: multiple center m/z values.
    * `Tuple{<:Real, <:Real}`: an m/z range.
    * `Nothing`: no m/z value is specified.
* `accuracy`: m/z accuracy, specified as a number or a [`Criteria`](@ref).
* `fwhm::Real`: full width at half maximum in m/z.
* `digits`: number of decimal places retained for m/z values. `nothing` means there is no limit.
* `stepsize`: step size for the mass scan. It must be a multiple of `10 ^ (-digits)`. If `digits` is `nothing`, `stepsize` must also be `nothing`.

# Constructors
    Quadrupole(mz = nothing; accuracy = 0.1, offset = 0, fwhm = 0.7, digits = nothing, stepsize = nothing, unit = 0.7, flatness = 100, taperproportion = 0.2)

## Arguments
* `offset::Real`: offset added to the window center; the resulting center is `mz + offset`.
* `unit::Real`: unit-mass resolution (FWHM), the largest FWHM that uses a Gaussian window. For larger FWHM values, a flat-topped super-Gaussian window is used.
* `flatness::Real`: flatness of the super-Gaussian window's Gaussian falloff. Practical values range from 10 to 1000.
* `taperproportion::Real`: proportion of the window occupied by the tapered region of a [`TukeyWindow`](@ref), used for single-m/z isolation. Values between 0 and 1 select a Tukey window; values of 1 or greater select a cosine window, and nonpositive values select a rectangular window.
"""
struct Quadrupole{W, M} <: AbstractMSAnalyzer{W, M}
    window::W
    mz::M
    accuracy
    fwhm
    digits
    stepsize
end 

function Quadrupole(mz = nothing; accuracy = 0.1, offset = 0, fwhm = 0.7, digits = nothing, stepsize = nothing, unit = 0.7, flatness = 100, taperproportion = 0.2)
    mz = isnothing(mz) ? mz : @.(mz + offset)
    if isnothing(mz) || mz isa Tuple
        if fwhm > unit
            f = sqrt(log(2, flatness))
            n = log(2, f) / log(2, unit / fwhm * (f - 1) + 1) * 2
            Quadrupole(SuperGaussianWindow(n), mz, accuracy, fwhm, digits, stepsize)
        else
            Quadrupole(GaussianWindow(), mz, accuracy, fwhm, digits, stepsize)
        end
    elseif taperproportion >= 1
        Quadrupole(CosineWindow(), mz, accuracy, fwhm, digits, stepsize)
    elseif taperproportion > 0
        Quadrupole(TukeyWindow(taperproportion), mz, accuracy, fwhm, digits, stepsize)
    else
        Quadrupole(RectWindow(), mz, accuracy, fwhm, digits, stepsize)
    end 
end

doc_qit = """
    QuadrupoleIonTrap{W, M} <: AbstractMSAnalyzer{W, M} 
    const QIT = QuadrupoleIonTrap

Quadrupole ion trap (QIT) mass analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value or values:
    * `Real`: one center m/z value.
    * `Vector{<:Real}`: multiple center m/z values.
    * `Tuple{<:Real, <:Real}`: an m/z range.
    * `Nothing`: no m/z value is specified.
* `accuracy`: m/z accuracy, specified as a number or a [`Criteria`](@ref).
* `fwhm::Real`: full width at half maximum in m/z.
* `digits`: number of decimal places retained for m/z values. `nothing` means there is no limit.
* `stepsize`: step size for the mass scan. It must be a multiple of `10 ^ (-digits)`. If `digits` is `nothing`, `stepsize` must also be `nothing`.

# Constructors
    QuadrupoleIonTrap(mz = nothing; accuracy = 0.1, offset = 0, fwhm = 0.7, stepsize = nothing, power = 4)

## Arguments
* `offset::Real`: offset added to the window center; the resulting center is `mz + offset`.
* `power::Real`: exponent of the [`SuperGaussianWindow`](@ref) used for SWIFT isolation. Practical values range from 4 to 8.
"""
@doc doc_qit
struct QuadrupoleIonTrap{W, M} <: AbstractMSAnalyzer{W, M} 
    window::W
    mz::M
    accuracy
    fwhm
    digits
    stepsize
end 

function QuadrupoleIonTrap(mz = nothing; accuracy = 0.1, offset = 0, fwhm = 0.7, stepsize = nothing, power = 4)
    mz = isnothing(mz) ? mz : @.(mz + offset)
    if isnothing(mz) || mz isa Tuple
        # Instability scan mode
        QuadrupoleIonTrap(GaussianWindow(), mz, accuracy, fwhm, digits, stepsize)
    else
        # SWIFT
        QuadrupoleIonTrap(SuperGaussianWindow(power), mz, accuracy, fwhm, digits, stepsize)
    end
end
@doc doc_qit
const QIT = QuadrupoleIonTrap

doc_lit = """
    LinearIonTrap{W, M} <: AbstractMSAnalyzer{W, M} 
    const LIT = LinearIonTrap

Linear ion trap (LIT) mass analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value or values:
    * `Real`: one center m/z value.
    * `Vector{<:Real}`: multiple center m/z values.
    * `Tuple{<:Real, <:Real}`: an m/z range.
    * `Nothing`: no m/z value is specified.
* `accuracy`: m/z accuracy, specified as a number or a [`Criteria`](@ref).
* `fwhm::Real`: full width at half maximum in m/z.
* `digits`: number of decimal places retained for m/z values. `nothing` means there is no limit.
* `stepsize`: step size for the mass scan. It must be a multiple of `10 ^ (-digits)`. If `digits` is `nothing`, `stepsize` must also be `nothing`.

# Constructors
    LinearIonTrap(mz = nothing; accuracy = 0.1, offset = 0, fwhm = 0.7, stepsize = nothing, taperproportion = 0.2)

## Arguments
* `offset::Real`: offset added to the window center; the resulting center is `mz + offset`.
* `taperproportion::Real`: proportion of the window occupied by the tapered region of a [`TukeyWindow`](@ref), used for SWIFT isolation. Values between 0 and 1 select a Tukey window; values of 1 or greater select a cosine window, and nonpositive values select a rectangular window.
"""
@doc doc_lit
struct LinearIonTrap{W, M} <: AbstractMSAnalyzer{W, M} 
    window::W
    mz::M
    accuracy
    fwhm
    digits
    stepsize
end 

function LinearIonTrap(mz = nothing; accuracy = 0.1, offset = 0, fwhm = 0.7, stepsize = nothing, taperproportion = 0.2)
    mz = isnothing(mz) ? mz : @.(mz + offset)
    if isnothing(mz) || mz isa Tuple
        # Instability scan mode
        LinearIonTrap(GaussianWindow(), mz, accuracy, fwhm, digits, stepsize)
    elseif taperproportion >= 1
        LinearIonTrap(CosineWindow(), mz, accuracy, fwhm, digits, stepsize)
    elseif taperproportion > 0
        LinearIonTrap(TukeyWindow(taperproportion), mz, accuracy, fwhm, digits, stepsize)
    else
        LinearIonTrap(RectWindow(), mz, accuracy, fwhm, digits, stepsize)
    end 
end
@doc doc_lit
const LIT = LinearIonTrap

doc_tof = """
    TimeOfFlight{W, M} <: AbstractMSAnalyzer{W, M} 
    const TOF = TimeofFlight

Time-of-flight (TOF) mass analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value or values:
    * `Real`: one center m/z value.
    * `Vector{<:Real}`: multiple center m/z values.
    * `Tuple{<:Real, <:Real}`: an m/z range.
    * `Nothing`: no m/z value is specified.
* `accuracy`: m/z accuracy, specified as a number or a [`Criteria`](@ref).
* `resolution`: theoretical maximum resolving power.
* `mz50`: m/z value at which the resolving power is half the specified `resolution`.

# Constructors
    TimeOfFlight(mz = nothing; window = GaussianWindow(), accuracy = rcrit(1e-6), resolution = 80000, mz50 = 120)
"""
@doc doc_tof
struct TimeOfFlight{W, M} <: AbstractMSAnalyzer{W, M} 
    window::W
    mz::M
    accuracy
    resolution
    mz50
end 
TimeOfFlight(mz = nothing; window = GaussianWindow(), accuracy = rcrit(1e-6), resolution = 80000, mz50 = 120) = TimeOfFlight(window, mz, accuracy, resolution, mz50)
# 100000 .* sqrt.(x) ./ (8 .+ sqrt.(x)) .- 20000
@doc doc_tof
const TOF = TimeOfFlight

"""
    Orbitrap{W, M} <: AbstractMSAnalyzer{W, M} 

Orbitrap mass analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value or values:
    * `Real`: one center m/z value.
    * `Vector{<:Real}`: multiple center m/z values.
    * `Tuple{<:Real, <:Real}`: an m/z range.
    * `Nothing`: no m/z value is specified.
* `accuracy`: m/z accuracy, specified as a number or a [`Criteria`](@ref).
* `resolution`: resolving power at m/z 200.
* `mz50`: m/z value at which the resolving power is half the specified `resolution`.

# Constructors
    Orbitrap(mz = nothing; window = GaussianWindow(), accuracy = rcrit(1e-6), resolution = 480000, mz50 = 800)
"""
struct Orbitrap{W, M} <: AbstractMSAnalyzer{W, M} 
    window::W
    mz::M
    accuracy
    resolution
    mz50
end 
Orbitrap(mz = nothing; window = GaussianWindow(), accuracy = rcrit(1e-6), resolution = 480000, mz50 = 800) = Orbitrap(window, mz, accuracy, resolution, mz50)

doc_fticr = """
    FourierTransformIonCyclotronResonance{W, M} <: AbstractMSAnalyzer{W, M} 
    const FTICR = FourierTransformIonCyclotronResonance

Fourier-transform ion cyclotron resonance (FTICR) mass analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value or values:
    * `Real`: one center m/z value.
    * `Vector{<:Real}`: multiple center m/z values.
    * `Tuple{<:Real, <:Real}`: an m/z range.
    * `Nothing`: no m/z value is specified.
* `accuracy`: m/z accuracy, specified as a number or a [`Criteria`](@ref).
* `resolution`: resolving power at m/z 200.
* `mz50`: m/z value at which the resolving power is half the specified `resolution`.

# Constructors
    FourierTransformIonCyclotronResonance(mz = nothing; window = GaussianWindow(), accuracy = rcrit(2e-7), resolution = 20000000, mz50 = 400)
"""
@doc doc_fticr
struct FourierTransformIonCyclotronResonance{W, M} <: AbstractMSAnalyzer{W, M} 
    window::W
    mz::M
    accuracy
    resolution
    mz50
end 
FourierTransformIonCyclotronResonance(mz = nothing; window = GaussianWindow(), accuracy = rcrit(2e-7), resolution = 20000000, mz50 = 400) = FourierTransformIonCyclotronResonance(window, mz, accuracy, resolution, mz50)
@doc doc_fticr
const FTICR = FourierTransformIonCyclotronResonance

abstract type AlgSpecPeak end
struct LocalMaxima <: AlgSpecPeak
    threshold::Real
end 
LocalMaxima() = LocalMaxima(0.1)

"""
    CoelutingIsobars

A type containing settings and results for analyzing co-eluting isobars.

# Fields
* `elution::Vector{<:Pair}`: pairs of elution functions (such as `retentiontime`) and criteria.
* `msanalyzer::Vector{<:Pair}`: pairs of mass analyzers and abundance criteria.
* `target::Table`: target chemicals represented as [`Isotopomers`](@ref).
* `isobar::Table`: potential isobars represented as [`Isotopomers`](@ref).
* `tables::Vector`: isotopologue tables for the target chemicals.

Each entry in `msanalyzer` corresponds to one MS stage of the transition.
Stages whose corresponding analyzer is `nothing` are not filtered.
If `msanalyzer` has fewer entries than the transition has stages, stages without a corresponding analyzer are not filtered.

# Constructors
    CoelutingIsobars(elution, msanalyzer, table::Table; nfilter = 2, kwargs...)
    CoelutingIsobars(elution, msanalyzer, target::Table, isobar::Table; nfilter = 2)

## Arguments
* `table::Table`: table of target chemicals not represented as [`Isotopomers`](@ref). The `target` and `isobar` tables are generated from this table.
* `nfilter::Int`: number of filtering stages to apply. `0` disables filtering and leaves `tables` empty; `1` filters by `elution` only; `2` filters by both `elution` and `msanalyzer`.
* Other keyword arguments are passed to [`Isotopologues`](@ref) to generate `isobar` from `table`.
"""
struct CoelutingIsobars 
    elution::Vector{<:Pair}
    msanalyzer::Vector{<:Pair}
    target::Table 
    isobar::Table
    tables::Vector
end

function CoelutingIsobars(elution, msanalyzer, table::Table; nfilter = 2, kwargs...) 
    isobar = TandemIsotopologues(table; kwargs...)
    target = isobar[[findfirst(x -> ischemicalequal(x, y), isobar.Chemical) for y in table.Chemical]]
    CoelutingIsobars(elution, msanalyzer, target, isobar; nfilter)
end

function CoelutingIsobars(elution, msanalyzer, target::Table, isobar::Table; nfilter = 2) 
    ci = CoelutingIsobars(elution, msanalyzer, target, isobar, Any[])
    nfilter < 1 && return ci 
    elution_filter!(ci)
    nfilter -= 1
    nfilter < 1 && return ci 
    ms_filter!(ci)
end