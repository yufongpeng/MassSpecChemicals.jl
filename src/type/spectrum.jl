"""
    Spectrum{T}

Spectrum type.

# Fields
* `spectrum::Vector{T}`: spectrum intensity vector.
* `initial_mass::T`: the first m/z value of `spectrum`.
* `binsize::T`: the size of m/z bin, i.e. m/z difference of consecutive elements.
* `stepsize::Int`: number of bins for each data point of the actual spectrum. 
* `table::Table`: the source table of ions. 
"""
struct Spectrum{T}
    spectrum::Vector{T}
    initial_mass::T 
    binsize::T
    stepsize::Int
    table::Table
end

abstract type AbstractWindow end
abstract type ZeroCenteredtWindow end
struct GaussianWindow <: AbstractWindow end
struct GaussianDilatedWindow <: AbstractWindow
    taperwidth
    taperheight
    min_fwhm
end
GaussianDilatedWindow(taperwidth, taperheight) = GaussianDilatedWindow(taperwidth, taperheight, taperwidth / (1 - sqrt(log(2, 1 / taperheight))))
struct GaussianTailedUniformWindow <: AbstractWindow
    taperwidth
end
struct FixedTaperTukeyWindow <: AbstractWindow 
    taperwidth
end
struct TukeyWindow <: AbstractWindow 
    taperproportion
end
struct CosineWindow <: AbstractWindow end
struct RectWindow <: AbstractWindow end
struct DiscreteGaussianWindow <: ZeroCenteredtWindow end
struct SampledGaussianWindow <: ZeroCenteredtWindow end
struct SuperGaussianWindow{T} <: AbstractWindow 
    power::T
end
struct PowerCosineWindow{T} <: AbstractWindow 
    power::T
end

abstract type AbstractMSAnalyzer{W, M} end 

"""
    MSAnalyzer{W, M} <: AbstractMSAnalyzer{W, M} 

Generic MS Analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value(s).
* `accuracy`: m/z accuracy, either a number or a `Criteria`.
* `fwhm_mz`: a function mapping m/z value to fwhm.
* `digits`: digits of m/z value.
* `stepsize`: step size for discrete mass scan such as quadrupole. It must be multiples of `10 ^ (-digits)`.
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

Quadrupole MS Analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value(s).
* `accuracy`: m/z accuracy, either a number or a `Criteria`.
* `fwhm::Real`: m/z fwhm.
* `digits`: digits of m/z value. 
* `stepsize`: step size for discrete mass scan such as quadrupole. It must be multiples of `10 ^ (-digits)`.

# Constructors
`Quadrupole([mz = nothing]; scan = true, accuracy = 0.1, offset = 0, fwhm = 0.7, digits = nothing, stepsize = nothing, unit = 0.7, flatness = 100, taperproportion = 0.2)`

## Arguments
* `scan`: whether operating in scan or SIM mode.
* `offset`: m/z offset of the center of window. The resulting m/z center will be `mz + offset`.
* `unit`: unit mass resolution (fwhm); the largest fwhm that remains gaussian. Window function becomes super gaussian (flat top) for fwhm larger than unit.
* `flatness`: how flat is the gaussian fall-off of the super gaussian function. Practical values range from 10 to 1000. 
* `taperproportion`: the proportion of the window that is tapered for `TukeyWindow` (used in single m/z isolation). Practical values range from 0.1 to 8. 
"""
struct Quadrupole{W, M} <: AbstractMSAnalyzer{W, M}
    window::W
    mz::M
    accuracy
    fwhm
    digits
    stepsize
end 

function Quadrupole(mz = nothing; scan = true, accuracy = 0.1, offset = 0, fwhm = 0.7, digits = nothing, stepsize = nothing, unit = 0.7, flatness = 100, taperproportion = 0.2)
    mz = isnothing(mz) ? mz : @.(mz + offset)
    if (isnothing(mz) || mz isa Tuple) && scan
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

"""
    QuadrupoleIonTrap{W, M} <: AbstractMSAnalyzer{W, M} 
    const QIT = QuadrupoleIonTrap

Quadrupole Ion Trap (QIT) MS Analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value(s).
* `accuracy`: m/z accuracy, either a number or a `Criteria`.
* `fwhm::Real`: m/z fwhm.
* `digits`: digits of m/z value. 
* `stepsize`: step size for discrete mass scan such as quadrupole. It must be multiples of `10 ^ (-digits)`.

# Constructors
`QuadrupoleIonTrap([mz = nothing]; scan = true, accuracy = 0.1, offset = 0, fwhm = 0.7, stepsize = nothing, power = 4)`

## Arguments
* `scan`: whether operating in scan or SIM mode.
* `offset`: m/z offset of the center of window. The resulting m/z center will be `mz + offset`.
* `power`: the power of super gaussian function (used in SWIFT isolation). Practical values range from 4 to 8. 
"""
struct QuadrupoleIonTrap{W, M} <: AbstractMSAnalyzer{W, M} 
    window::W
    mz::M
    accuracy
    fwhm
    digits
    stepsize
end 

function QuadrupoleIonTrap(mz = nothing; scan = true, accuracy = 0.1, offset = 0, fwhm = 0.7, stepsize = nothing, power = 4)
    mz = isnothing(mz) ? mz : @.(mz + offset)
    if (isnothing(mz) || mz isa Tuple) && scan
        # Instability scan mode
        QuadrupoleIonTrap(GaussianWindow(), mz, accuracy, fwhm, digits, stepsize)
    else
        # SWIFT
        QuadrupoleIonTrap(SuperGaussianWindow(power), mz, accuracy, fwhm, digits, stepsize)
    end
end

const QIT = QuadrupoleIonTrap

"""
    LinearIonTrap{W, M} <: AbstractMSAnalyzer{W, M} 
    const LIT = LinearIonTrap

Linear Ion Trap (LIT) MS Analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value(s).
* `accuracy`: m/z accuracy, either a number or a `Criteria`.
* `fwhm::Real`: m/z fwhm.
* `digits`: digits of m/z value. 
* `stepsize`: step size for discrete mass scan such as quadrupole. It must be multiples of `10 ^ (-digits)`.

# Constructors
`LinearIonTrap([mz = nothing]; scan = true, accuracy = 0.1, offset = 0, fwhm = 0.7, stepsize = nothing, taperproportion = 0.2)`

## Arguments
* `scan`: whether operating in scan or SIM mode.
* `offset`: m/z offset of the center of window. The resulting m/z center will be `mz + offset`.
* `taperproportion`: the proportion of the window that is tapered for `TukeyWindow` (used in SWIFT isolation).
"""
struct LinearIonTrap{W, M} <: AbstractMSAnalyzer{W, M} 
    window::W
    mz::M
    accuracy
    fwhm
    digits
    stepsize
end 

function LinearIonTrap(mz = nothing; scan = true, accuracy = 0.1, offset = 0, fwhm = 0.7, stepsize = nothing, taperproportion = 0.2)
    mz = isnothing(mz) ? mz : @.(mz + offset)
    if (isnothing(mz) || mz isa Tuple) && scan 
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

const LIT = LinearIonTrap

"""
    TimeOfFlight{W, M} <: AbstractMSAnalyzer{W, M} 
    const TOF = TimeofFlight

Time-of-Flight (TOF) MS Analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value(s).
* `accuracy`: m/z accuracy, either a number or a `Criteria`.
* `resolution`: theoretical maximal m/z resolution. 
* `mz50`: m/z value with half `resolution`.

# Constructors
`TimeOfFlight([mz = nothing]; window = GaussianWindow(), accuracy = rcrit(1e-6), resolution = 66000, mz50 = 120)`
"""
struct TimeOfFlight{W, M} <: AbstractMSAnalyzer{W, M} 
    window::W
    mz::M
    accuracy
    resolution
    mz50
end 
TimeOfFlight(mz = nothing; window = GaussianWindow(), accuracy = rcrit(1e-6), resolution = 66000, mz50 = 120) = TimeOfFlight(window, mz, accuracy, resolution, mz50)
# 100000 .* sqrt.(x) ./ (8 .+ sqrt.(x)) .- 20000
const TOF = TimeOfFlight

"""
    Orbitrap{W, M} <: AbstractMSAnalyzer{W, M} 

Orbitrap MS Analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value(s).
* `accuracy`: m/z accuracy, either a number or a `Criteria`.
* `resolution`: m/z resolution @200. 
* `mz50`: m/z value with half `resolution`.

# Constructors
`Orbitrap([mz = nothing]; window = GaussianWindow(), accuracy = rcrit(1e-6), resolution = 480000, mz50 = 800)`
"""
struct Orbitrap{W, M} <: AbstractMSAnalyzer{W, M} 
    window::W
    mz::M
    accuracy
    resolution
    mz50
end 
Orbitrap(mz = nothing; window = GaussianWindow(), accuracy = rcrit(1e-6), resolution = 480000, mz50 = 800) = Orbitrap(window, mz, accuracy, resolution, mz50)

"""
    FourierTransformIonCyclotronResonance{W, M} <: AbstractMSAnalyzer{W, M} 
    const FTICR = FourierTransformIonCyclotronResonance

Fourier-transform ion cyclotron resonance (FTICR) MS Analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value(s).
* `accuracy`: m/z accuracy, either a number or a `Criteria`.
* `resolution`: m/z resolution @200. 
* `mz50`: m/z value with half `resolution`.

# Constructors
`FourierTransformIonCyclotronResonance(mz = nothing; window = GaussianWindow(), accuracy = rcrit(2e-7), resolution = 20000000, mz50 = 400)`
"""
struct FourierTransformIonCyclotronResonance{W, M} <: AbstractMSAnalyzer{W, M} 
    window::W
    mz::M
    accuracy
    resolution
    mz50
end 
FourierTransformIonCyclotronResonance(mz = nothing; window = GaussianWindow(), accuracy = rcrit(2e-7), resolution = 20000000, mz50 = 400) = FourierTransformIonCyclotronResonance(window, mz, accuracy, resolution, mz50)
const FTICR = FourierTransformIonCyclotronResonance

abstract type AlgSpecPeak end
struct LocalMaxima <: AlgSpecPeak end 
struct FWHMMaxima <: AlgSpecPeak end 