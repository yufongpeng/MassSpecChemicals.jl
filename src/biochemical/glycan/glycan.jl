module Glycans

using Reexport
@reexport using ..BioChemicals
using ..MSChemicals: AbstractChemical
export Saccharide, Monosaccharide, 
        AbstractHex, Hex, Glc, Gal, Man, 
        AbstractdHex, dHex, Fuc, 
        AbstractHexN, HexN, GlcN, GalN, ManN, 
        AbstractHexNAc, HexNAc, GlcNAc, GalNAc, ManNAc, 
        AbstractHexA, HexA, GlcA, GalA, ManA,
        SialicAcid, Neu, Neu5Ac, NeuAc, Neu5Gc, NeuGc, Kdb,
        AbstractPen, Pen, Rib, Ara, Xyl, AbstractdPen, dPen, dRib, 
        Inositol, Sulfoquinovose, Glycan, GlyComp,
        AbstractGlycan,
        AbstractAnomerposition, Anomerposition, Alphaposition, Betaposition,
        Ganglioseries,
        GM4,
        SM4,
        Lac,
        SM3,
        SM2,
        SM1a,
        SM1b,
        SB1a,
        Ganglio0series,
        GA2,
        GA1,
        GM1b,
        GM1α,
        GD1c,
        GD1α,
        GD1e,
        GanglioAseries,
        GM3,
        GM2,
        GM1a,
        GD1a,
        GD1aα,
        GT1a,
        GT1aα,
        GanglioBseries,
        GD3,
        GD2,
        GD1b,
        GT1b,
        GT1bα,
        GQ1b,
        GQ1bα,
        GanglioCseries,
        GT3,
        GT2,
        GT1c,
        GQ1c,
        GQ1cα,
        GP1c,
        GP1cα,
        GM1,
        GD1,
        GT1,
        GQ1,
        GP1,
        Globoseries,
        Gb3,
        Gb4,
        Isogloboseries,
        iGb3,
        iGb4,
        Lactoseries,
        Lc3,
        LM1
        
abstract type Saccharide <: AbstractChemical end
abstract type Monosaccharide{T} <: Saccharide end
abstract type AbstractHex{T} <: Monosaccharide{T} end
abstract type AbstractdHex{T} <: Monosaccharide{T} end
abstract type AbstractHexN{T} <: Monosaccharide{T} end
abstract type AbstractHexNAc{T} <: Monosaccharide{T} end
abstract type AbstractHexA{T} <: Monosaccharide{T} end
abstract type SialicAcid{T} <: Monosaccharide{T} end
abstract type AbstractPen{T} <: Monosaccharide{T} end
abstract type AbstractdPen{T} <: Monosaccharide{T} end
#=
T
Nothing
Vector{<: Union{FunctionalGroup, AbstractChemical}}
Vector{<: Pair{UInt8, <: Union{FunctionalGroup, AbstractChemical}}
=#
"""
    Hex{T} <: AbstractHex{T}

Hexose with or without substituents.
"""
struct Hex{T} <: AbstractHex{T}
    substituent::T
end
Hex() = Hex(nothing)
"""
    Glc{T} <: AbstractHex{T}

Glucose with or without substituents.
"""
struct Glc{T} <: AbstractHex{T}
    substituent::T
end
Glc() = Glc(nothing)
"""
    Gal{T} <: AbstractHex{T}

Galactose with or without substituents.
"""
struct Gal{T} <: AbstractHex{T}    
    substituent::T
end
Gal() = Gal(nothing)
"""
    Man{T} <: AbstractHex{T}

Mannose with or without substituents.
"""
struct Man{T} <: AbstractHex{T}    
    substituent::T
end
Man() = Man(nothing)

"""
    HexN{T} <: AbstractHexN{T}

Hexosamine with or without substituents.
"""
struct HexN{T} <: AbstractHexN{T}
    substituent::T
end
HexN() = HexN(nothing)
"""
    GlcN{T} <: AbstractHexN{T}

Glucosamine with or without substituents.
"""
struct GlcN{T} <: AbstractHexN{T}
    substituent::T
end
GlcN() = GlcN(nothing)
"""
    GalN{T} <: AbstractHexN{T}

Galactosamine with or without substituents.
"""
struct GalN{T} <: AbstractHexN{T}    
    substituent::T
end
GalN() = GalN(nothing)
"""
    ManN{T} <: AbstractHexN{T}

Mannosamine with or without substituents.
"""
struct ManN{T} <: AbstractHexN{T}    
    substituent::T
end
ManN() = ManN(nothing)

"""
    HexNAc{T} <: AbstractHexNAc{T}

N-acetylhexosamine with or without substituents.
"""
struct HexNAc{T} <: AbstractHexNAc{T}
    substituent::T
end
HexNAc() = HexNAc(nothing)
"""
    GlcNAc{T} <: AbstractHexNAc{T}

N-acetylglucosamine with or without substituents.
"""
struct GlcNAc{T} <: AbstractHexNAc{T}
    substituent::T
end
GlcNAc() = GlcNAc(nothing)
"""
    GalNAc{T} <: AbstractHexNAc{T}

N-acetylgalactosamine with or without substituents.
"""
struct GalNAc{T} <: AbstractHexNAc{T}    
    substituent::T
end
GalNAc() = GalNAc(nothing)
"""
    ManNAc{T} <: AbstractHexNAc{T}

N-acetylmannosamine with or without substituents.
"""
struct ManNAc{T} <: AbstractHexNAc{T}    
    substituent::T
end
ManNAc() = ManNAc(nothing)

"""
    HexA{T} <: AbstractHexA{T}

Hexuronate with or without substituents.
"""
struct HexA{T} <: AbstractHexA{T}
    substituent::T
end
HexA() = HexA(nothing)
"""
    GlcA{T} <: AbstractHexA{T}

Glucuronic acid with or without substituents.
"""
struct GlcA{T} <: AbstractHexA{T}
    substituent::T
end
GlcA() = GlcA(nothing)
"""
    GalA{T} <: AbstractHexA{T}

Galactcuronic acid with or without substituents.
"""
struct GalA{T} <: AbstractHexA{T}
    substituent::T
end
GalA() = GalA(nothing)
"""
    ManA{T} <: AbstractHexA{T}

Mannuronic acid with or without substituents.
"""
struct ManA{T} <: AbstractHexA{T}
    substituent::T
end
ManA() = ManA(nothing)

"""
    dHex{T} <: AbstractdHex{T}

Deoxyhexose with or without substituents.
"""
struct dHex{T} <: AbstractdHex{T}
    substituent::T
end
dHex() = dHex(nothing)
"""
    Fuc{T} <: AbstractdHex{T}

Fucose with or without substituents.
"""
struct Fuc{T} <: AbstractdHex{T}
    substituent::T
end
Fuc() = Fuc(nothing)

"""
    Neu{T} <: SialicAcid{T}

Neuraminic acid with or without substituents.
"""
struct Neu{T} <: SialicAcid{T}
    substituent::T
end
Neu() = Neu(nothing)
"""
    Neu5Ac{T} <: SialicAcid{T}

N-acetylneuraminic acid (nana) with or without substituents.
"""
struct Neu5Ac{T} <: SialicAcid{T}
    substituent::T
end
const NeuAc = Neu5Ac
Neu5Ac() = Neu5Ac(nothing)
"""
    Neu5Gc{T} <: SialicAcid{T}

N-glycolylneuraminic acid with or without substituents.
"""
struct Neu5Gc{T} <: SialicAcid{T}
    substituent::T
end
const NeuGc = Neu5Gc
Neu5Gc() = Neu5Gc(nothing)
"""
    Kdn <: SialicAcid

Keto-deoxy-glycero-glactonononic acid with or without substituents.
"""
struct Kdb{T} <: SialicAcid{T}
    substituent::T
end
Kdb() = Kdb(nothing)

"""
    Pen{T} <: AbstractPen{T}

Pentose with or without substituents.
"""
struct Pen{T} <: AbstractPen{T}
    substituent::T
end
Pen() = Pen(nothing)
"""
    Ara{T} <: AbstractPen{T}

Arabinose with or without substituents.
"""
struct Ara{T} <: AbstractPen{T}
    substituent::T
end
Ara() = Ara(nothing)
"""
    Xyl{T} <: AbstractPen{T}

Xylose with or without substituents.
"""
struct Xyl{T} <: AbstractPen{T}
    substituent::T
end
Xyl() = Xyl(nothing)
"""
    Rib{T} <: AbstractPen{T}

Ribose with or without substituents.
"""
struct Rib{T} <: AbstractPen{T}
    substituent::T
end
Rib() = Rib(nothing)

"""
    dPen{T} <: AbstractdPen{T}

Deoxypentose with or without substituents.
"""
struct dPen{T} <: AbstractdPen{T}
    substituent::T
end
dPen() = dPen(nothing)
"""
    dRib{T} <: AbstractdPen{T}

Deoxyribose with or without substituents.
"""
struct dRib{T} <: AbstractdPen{T}
    substituent::T
end
dRib() = dRib(nothing)

struct Inositol{T} <: Monosaccharide{T}
    substituent::T
end
Inositol() = Inositol(nothing)

struct Sulfoquinovose{T} <: Monosaccharide{T}
    substituent::T
end
Sulfoquinovose() = Sulfoquinovose(nothing)

abstract type AbstractGlycan <: Saccharide end
struct GM4 <: AbstractGlycan end
struct SM4 <: AbstractGlycan end
struct Lac <: AbstractGlycan end
abstract type Ganglioseries <: AbstractGlycan end
struct SM3 <: Ganglioseries end
struct SM2 <: Ganglioseries end
struct SM1a <: Ganglioseries end
struct SM1b <: Ganglioseries end
struct SB1a <: Ganglioseries end
abstract type Ganglio0series <: Ganglioseries end
struct GA2 <: Ganglio0series end
struct GA1 <: Ganglio0series end
struct GM1b <: Ganglio0series end
struct GM1α <: Ganglio0series end
struct GD1c <: Ganglio0series end
struct GD1α <: Ganglio0series end
const GD1e = GD1α
abstract type GanglioAseries <: Ganglioseries end
struct GM3 <: GanglioAseries end
struct GM2 <: GanglioAseries end
struct GM1a <: GanglioAseries end
struct GD1a <: GanglioAseries end
struct GD1aα <: GanglioAseries end
struct GT1a <: GanglioAseries end
struct GT1aα <: GanglioAseries end
abstract type GanglioBseries <: Ganglioseries end
struct GD3 <: GanglioBseries end
struct GD2 <: GanglioBseries end
struct GD1b <: GanglioBseries end
struct GT1b <: GanglioBseries end
struct GT1bα <: GanglioBseries end
struct GQ1b <: GanglioBseries end
struct GQ1bα <: GanglioBseries end
abstract type GanglioCseries <: Ganglioseries end
struct GT3 <: GanglioCseries end
struct GT2 <: GanglioCseries end
struct GT1c <: GanglioCseries end
struct GQ1c <: GanglioCseries end
struct GQ1cα <: GanglioCseries end
struct GP1c <: GanglioCseries end
struct GP1cα <: GanglioCseries end
struct GM1 <: Ganglioseries
    code::UInt8
end # GM1a, x, x, GM1b, x, x, x, GM1α
struct GD1 <: Ganglioseries
    code::UInt8
end # GD1a, GD1b, x, GD1c, GD1aα, x, x, GD1α 
struct GT1 <: Ganglioseries
    code::UInt8
end # GT1a, GT1b, GT1c, x, GT1aα, GT1bα, x, x
struct GQ1 <: Ganglioseries
    code::UInt8
end # x, GQ1b, GQ1c, x, x, GQ1bα, GQ1cα, x
struct GP1 <: Ganglioseries
    code::UInt8
end # x, x GP1c, x, x, x, GP1cα, x
abstract type Globoseries <: AbstractGlycan end
struct Gb3 <: Globoseries end
struct Gb4 <: Globoseries end
abstract type Isogloboseries <: AbstractGlycan end
struct iGb3 <: Isogloboseries end
struct iGb4 <: Isogloboseries end
abstract type Lactoseries <: AbstractGlycan end
struct Lc3 <: Lactoseries end
struct LM1 <: Lactoseries end
abstract type AbstractAnomerposition <: AbstractLinkageposition end
struct Anomerposition <: AbstractAnomerposition
    position::UInt8
end
ap(x) = Anomerposition(UInt8(x))
struct Alphaposition <: AbstractAnomerposition
    position::UInt8
end
α(x) = Alphaposition(UInt8(x))
struct Betaposition <: AbstractAnomerposition
    position::UInt8
end
β(x) = Betaposition(UInt8(x))
"""
    Glycan{T} <: Saccharide

Oligosaccarides or Polysaccharides with or without defined order.
"""
struct Glycan{T} <: AbstractGlycan
    chain::T
    linkage::Vector{Pair{AbstractAnomerposition, UInt8}}
end
Glycan(sugar...) = Glycan(sugar, Pair{AbstractAnomerposition, UInt8}[])
#= Tuple: order, length of linkage == chain, defined reducing end anomer, 
S::UInt8

S::Pair{UInt8, UInt8} x, y => p, a = divrem(x, 3) ap-y, ay
a = 0, unknown anomer; a = 1, α; a = 2, β
p = 0, default reducing end 1/2
=#
push!(AbstractChainedChemical, Glycan)

struct GlyComp{T} <: AbstractGlycan
    comp::Vector{Pair{T, UInt8}}
end
# Vector{Pair{Monosaccharide, UInt8}}: type => number, no order
end