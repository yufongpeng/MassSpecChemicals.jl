module Lipids
using Reexport, IterTools, MLStyle
@reexport using ..BioChemicals
using ..MSChemicals: AbstractChemical
@reexport using ..BioChemicals.BasicCompounds, ..BioChemicals.Metabolites, ..BioChemicals.AminoAcids, ..BioChemicals.Glycans
import ..BioChemicals: originalmolecule, leavinggroup, conjugation
import ..MSChemicals: chemicalname
using ..BioChemicals: lk, makemolecule, makelinkage
using ..BioChemicals.Glycans: ap, α, β
export AbstractCarbonChain, CarbonChain, IsoprenoidChain, Acyl, Alkyl, Alkenyl, SPB, AbstractSTRing, STRing, SRing, DSMSRing, DCRing, CASRing, BRSRing, EGSRing, DEGSRing, SISRing, STSRing,
        Lipid, 
        FattyAcyl, Hydrocarbon, FattyAcid, FattyAldehyde, WaxEster, FattyAmide, NacylAlkylAmine, NacylAmine, FattyAcylCarnitine, FattyAcylCoA, FattyAcylEstolide,
        Glycerolipid, Monoradylglycerol, Diradylglycerol, Triradylglycerol, Estolide, 
        Omodifiedmonoradylglycerol, Sulfoquinovosylmonoradylglycerol, Monogalactosylmonoradylglycerol, Digalactosylmonoradylglycerol,
        Omodifieddiradylglycerol, Sulfoquinovosyldiradylglycerol, Monogalactosyldiradylglycerol, Digalactosyldiradylglycerol,

        Glycerophospholipid, GlycerophosphoNacylethanolamine, 
        Phosphatidicacid, Phosphatidylcholine, Phosphatidylethanolamine, PhosphatidylNmethylethanolamine, PhosphatidylNNdimethylethanolamine, Phosphatidylserine, Phosphatidylinositol, Phosphatidylglycerol, Phosphatidylmethanol, Phosphatidylethanol,
        AbstractPhosphatidylinositolphosphate, Phosphatidylinositolphosphate, Phosphatidylinositolbiphosphate, Phosphatidylinositoltriphosphate, Phosphatidylglycerolphosphate, PhosphatidylNmodifiedethanolamine, PhosphatidylNmodifiedserine, Phosphatidylglycan,
        Lysophosphatidicacid, Lysophosphatidylcholine, Lysophosphatidylethanolamine, LysophosphatidylNmethylethanolamine, LysophosphatidylNNdimethylethanolamine, Lysophosphatidylserine, Lysophosphatidylinositol, Lysophosphatidylglycerol, Lysophosphatidylmethanol, Lysophosphatidylethanol, 
        AbstractLysophosphatidylinositolphosphate, Lysophosphatidylinositolphosphate, Lysophosphatidylinositolbiphosphate, Lysophosphatidylinositoltriphosphate, Lysophosphatidylglycerolphosphate, LysophosphatidylNmodifiedethanolamine, LysophosphatidylNmodifiedserine, Lysophosphatidylglycan,
        Lysobisphosphatidicacid, Semilysobisphosphatidicacid, Bisphosphatidicacid, Dilysocardiolipin, Monolysocardiolipin, Cardiolipin,

        Sphingolipid, Ceramide, SphingoidBase, Glycosylceramide, Glycosylsphingoidbase,
        CeramidePhosphate, SphingoidBasePhosphate, Sphingomyelin, Lysosphingomyelin, Inositolphosphorylceramide, Ethanolaminephosphorylceramide, Glycosylinositolphosphorylceramide, Mannosylinositolphosphorylceramide, Mannosyldiinositolphosphorylceramide,
        Lysoinositolphosphorylceramide, Lysoethanolaminephosphorylceramide, Lysoglycosylinositolphosphorylceramide, Lysomannosylinositolphosphorylceramide, Lysomannosyldiinositolphosphorylceramide,
        Sulfonolipid, Lysosulfonolipid, Acylceramide, MixSphingoBone, Acylhexosylceramide, Acylsphingomyelin, 

        Sterol, FreeSterol, Sterylester, SubstitutedSterol,
        Prenol, Retinylester, CoenzymeQ

include("utils.jl")
abstract type AbstractCarbonChain <: FunctionalGroup{Nothing, Nothing} end
struct CarbonChain{T, D, S, I <: Union{Nothing, String}} <: AbstractCarbonChain
    carbon::UInt8
    doublebond::D
    substituent::S
    isotopiclabel::I
end
CarbonChain{T}(carbon::C, doublebond, substituent, isotopiclabel = nothing) where {C <: Number, T} = CarbonChain{T}(UInt8(carbon), uint8ize(doublebond), uint8ize(substituent), isotopiclabel)
function CarbonChain{T}(carbon::UInt8, 
        doublebond::D, 
        substituent::S, 
        isotopiclabel::I = nothing
    ) where {D <: Union{UInt8, Vector{UInt8}, Vector{<: Pair{UInt8, Bool}}}, S <: Union{Nothing, UInt8, Vector{<: Pair{<: FunctionalGroup, UInt8}}, Vector{<: Pair{UInt8, <: FunctionalGroup}}}, I, T}
    CarbonChain{T, D, S, I}(carbon, doublebond, substituent, isotopiclabel)
end

struct IsoprenoidChain{N, I <: Union{Nothing, String}} <: AbstractCarbonChain
    isotopiclabel::I
end

abstract type CarbonChainType end
abstract type Radyl <: CarbonChainType end
struct Acyl <: Radyl end
struct Alkyl <: Radyl end
struct Alkenyl <: Radyl end
struct SPB <: CarbonChainType end
abstract type AbstractSTRing <: CarbonChainType end
struct STRing <: AbstractSTRing end
struct CRing <: AbstractSTRing end # Cholesterol
struct DSMSRing <: AbstractSTRing end # Desmosterol
struct DCRing <: AbstractSTRing end # Dihydrocholesterol
struct CASRing <: AbstractSTRing end # Campesterol
struct BRSRing <: AbstractSTRing end # Brassicasterol
struct EGSRing <: AbstractSTRing end # Ergosterol
struct DEGSRing <: AbstractSTRing end # Dehydroergosterol
struct SISRing <: AbstractSTRing end # Sitosterol
struct STSRing <: AbstractSTRing end # Stigmasterol
struct BAing <: AbstractSTRing end # Bile acid

#= 
T
CarbonChainType
Tuple: order STRing, SPB, Alkenyl Alkyl Acyl, more -> few
D
UInt8: ndoublebond
Vector{UInt8}: doublebond position
Vector{Pair{UInt8, Bool}}: doublebond position => Z(false)/E(true)

F
UInt8: noxygen
Vector{Pair{FunctionalGroup, UInt8}}: FunctionalGroup => number
Vector{Pair{UInt8, FunctionalGroup}}: Position => FunctionalGroup

C 
CarbonChain, Tuple{<: CarbonChain}
=#
const AlkylAcylChain = Union{<: CarbonChain{Tuple{<: Union{Alkyl, Alkenyl}, <: Acyl}}, Tuple{<: CarbonChain{<: Union{Alkyl, Alkenyl}}, <: CarbonChain{Acyl}}}
const Acyl2Chain = Union{<: CarbonChain{<: Tuple{<: Acyl, <: Acyl}}, Tuple{<: CarbonChain{Acyl}, <: CarbonChain{Acyl}}}
const Alkyl2AcylChain = Union{<: CarbonChain{Tuple{<: Union{Alkyl, Alkenyl}, <: Union{Alkyl, Alkenyl}, <: Acyl}}, Tuple{<: CarbonChain{<: Tuple{<: Union{Alkyl, Alkenyl}, <: Union{Alkyl, Alkenyl}}}, <: CarbonChain{Acyl}}, Tuple{<: CarbonChain{<: Union{Alkyl, Alkenyl}}, <: CarbonChain{<: Tuple{<: Union{Alkyl, Alkenyl}, Acyl}}}, Tuple{<: CarbonChain{<: Union{Alkyl, Alkenyl}}, <: CarbonChain{<: Union{Alkyl, Alkenyl}}, <: CarbonChain{Acyl}}}
abstract type Lipid{B, C} <: AbstractChemical end
struct FattyAcyl{B, C <: CarbonChain} <: Lipid{B, C}
    backbone::B
    chain::C
end

const Hydrocarbon{B, C} = FattyAcyl{B, C} where {B <: includeSIL(Dihydrogen), C <: Union{<: CarbonChain{Alkyl}, <: CarbonChain{Alkenyl}}}
const FattyAcid{B, C} = FattyAcyl{B, C} where {B <: includeSIL(HydrogenOxide), C <: CarbonChain{Acyl}}
const FattyAldehyde{B, C} = FattyAcyl{B, C} where {B <: includeSIL(Dihydrogen), C <: CarbonChain{Acyl}} # FAL
const FattyAlcohol{B, C} = FattyAcyl{B, C} where {B <: includeSIL(HydrogenOxide), C <: Union{<: CarbonChain{Alkyl}, <: CarbonChain{Alkenyl}}} # FOH
const WaxEster{B, C} = FattyAcyl{B, C} where {B <: includeSIL(HydrogenOxide), C <: AlkylAcylChain}
const FattyAmide{B, C} = FattyAcyl{B, C} where {B <: includeSIL(Ammonia), C <: CarbonChain{Acyl}}
const NacylAlkylAmine{B, C} = FattyAcyl{B, C} where {B <: includeSIL(Ammonia), C <: Union{<: AlkylAcylChain, <: Alkyl2AcylChain}} # NA
const NacylAmine{B, C} = FattyAcyl{B, C} where {B, C <: CarbonChain{Acyl}} # NAE, NAGly, NASer, NAT...
const FattyAcylCarnitine{B, C} = FattyAcyl{B, C} where {B <: includeSIL(Carnitine), C <: CarbonChain{Acyl}}  
const FattyAcylCoA{B, C} = FattyAcyl{B, C} where {B <: includeSIL(CoA), C <: CarbonChain{Acyl}}
const FattyAcylEstolid{B, C} = FattyAcyl{B, C} where {B <: includeSIL(HydrogenOxide), C <: Acyl2Chain}

const MonoradylChain = CarbonChain{<: Radyl}
const DiradylChain = Union{<: CarbonChain{Tuple{<: Radyl, <: Radyl}}, <: Tuple{<: CarbonChain{<: Radyl}, <: CarbonChain{<: Radyl}}}
# divrem 3
const MonoDiradylChain = Union{<: MonoradylChain, <: DiradylChain}
const TriradylChain = Union{<: CarbonChain{Tuple{<: Radyl, <: Radyl, <: Radyl}}, <: Tuple{<: CarbonChain{<: Radyl}, <: CarbonChain{<: Tuple{<: Radyl, <: Radyl}}}, <: Tuple{<: CarbonChain{<: Radyl}, <: CarbonChain{<: Radyl}, <: CarbonChain{<: Radyl}}}
const TetraradylChain = Union{<: CarbonChain{Tuple{<: Radyl, <: Radyl, <: Radyl, <: Radyl}}, <: Tuple{<: CarbonChain{<: Tuple{<: CarbonChain{<: Radyl}, <: Radyl, <: Radyl, <: Radyl}}}, 
                            <: Tuple{<: CarbonChain{<: Tuple{<: Radyl, <: Radyl}}, <: CarbonChain{<: Tuple{<: Radyl, <: Radyl}}}, 
                            <: Tuple{<: CarbonChain{<: CarbonChain{<: Radyl}, <: CarbonChain{<: Radyl}, <: Tuple{<: Radyl, <: Radyl}}}, 
                            <: Tuple{<: CarbonChain{<: Radyl}, <: CarbonChain{<: Radyl}, <: CarbonChain{<: Radyl}, <: CarbonChain{<: Radyl}}}

abstract type Glycerolipid{B, C} <: Lipid{B, C} end
struct Monoradylglycerol{B <: includeSIL(Glycerol), C <: MonoradylChain} <: Glycerolipid{B, C}
    backbone::B
    chain::C
    sn::UInt8
end
struct Diradylglycerol{B <: includeSIL(Glycerol), C <: DiradylChain} <: Glycerolipid{B, C}
    backbone::B
    chain::C
    sn::UInt8
end # DG [0] [0, 0] [1, 0] [1, 2] [1, 1]
struct Triradylglycerol{B <: includeSIL(Glycerol), C <: TriradylChain} <: Glycerolipid{B, C}
    backbone::B
    chain::C
    sn::UInt8
end # TG [0] [0, 0] [0, 0, 0] [1/2, 0] [1, 0, 0] [1, 2, 1] # divrem 16 divrem 4
struct Estolide{B <: includeSIL(Glycerol), C <: TriradylChain} <: Glycerolipid{B, C}
    backbone::B
    chain::C
    sn::UInt8
end
struct Omodifiedradylglycerol{B <: DehydratedChemical, C <: MonoDiradylChain} <: Glycerolipid{B, C}
    backbone::B
    chain::C
    sn::UInt8
end # DehydratedChemical((Glycerol(), ?), [[Odehydrogen(), ?], [1, ?]])
const Omodifiedmonoradylglycerol{B, C} = Omodifiedradylglycerol{B, C} where {B, C <: MonoradylChain}
const Omodifieddiradylglycerol{B, C} = Omodifiedradylglycerol{B, C} where {B, C <: DiradylChain}
const Sulfoquinovosylradylglycerol{B, C} = Omodifiedradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(Glycerol), <: includeSIL(Glycan{<: Tuple{<: Sulfoquinovose}})}}, C}
const Sulfoquinovosylmonoradylglycerol{B, C} = Sulfoquinovosylradylglycerol{B, C} where {B, C <: MonoradylChain}
const Sulfoquinovosyldiradylglycerol{B, C} = Sulfoquinovosylradylglycerol{B, C} where {B, C <: DiradylChain}
const Monogalactosylradylglycerol{B, C} = Omodifiedradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(Glycerol), <: includeSIL(Glycan{<: Tuple{<: Gal}})}}, C}
const Monogalactosylmonoradylglycerol{B, C} = Monogalactosylradylglycerol{B, C} where {B, C <: MonoradylChain}
const Monogalactosyldiradylglycerol{B, C} = Monogalactosylradylglycerol{B, C} where {B, C <: DiradylChain}
const Digalactosylradylglycerol{B, C} = Omodifiedradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(Glycerol), <: includeSIL(Glycan{<: Tuple{<: Gal, <: Gal}})}}, C}
const Digalactosylmonoradylglycerol{B, C} = Digalactosylradylglycerol{B, C} where {B, C <: MonoradylChain}
const Digalactosyldiradylglycerol{B, C} = Digalactosylradylglycerol{B, C} where {B, C <: DiradylChain}
const Glucuronosylradylglycerol{B, C} = Omodifiedradylglycerol{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(Glycerol), <: includeSIL(Glycan{<: Tuple{<: GlcA}})}}, C}
const Glucuronosylmonoradylglycerol{B, C} = Glucuronosylradylglycerol{B, C} where {B, C <: MonoradylChain}
const Glucuronosyldiradylglycerol{B, C} = Glucuronosylradylglycerol{B, C} where {B, C <: DiradylChain}
# SMGDG, CDPDAG
abstract type Glycerophospholipid{B, C} <: Lipid{B, C} end # PX [0] [0, 0] [1, 2], LPX [0] [1/2]
struct Radylglycerophosphate{B <: DehydratedChemical, C} <: Glycerophospholipid{B, C}
    backbone::B
    chain::C
    sn::UInt8
end

const PPA = includeSIL(PhosphoricAcid)
glycerophospho(T) = DehydratedChemical{<: Tuple{<: includeSIL(Glycerol), <: PPA, <: Union{<: T, <: IsotopiclabeledChemical{<: T}}}}
glycerophospho(T, S) = DehydratedChemical{<: Tuple{<: includeSIL(Glycerol), <: PPA, <: Union{<: T, <: IsotopiclabeledChemical{<: T}}, <: Union{<: S, <: IsotopiclabeledChemical{<: S}}}}
glycerophosphophosphate(T) = DehydratedChemical{<: Tuple{<: includeSIL(Glycerol), <: PPA, <: Union{<: T, <: IsotopiclabeledChemical{<: T}}, <: PPA}}
const Monoradylglycerophosphate{B, C} = Radylglycerophosphate{B, C} where {B, C <: MonoradylChain}
const Diradylglycerophosphate{B, C} = Radylglycerophosphate{B, C} where {B, C <: DiradylChain}
const Lysophosphatidicacid{B, C} = Monoradylglycerophosphate{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(Glycerol), <: PPA}}, C}
const Phosphatidicacid{B, C} = Diradylglycerophosphate{B, C} where {B <: DehydratedChemical{<: Tuple{<: includeSIL(Glycerol), <: PPA}}, C}
const Lysophosphatidylcholine{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Choline), C}
const Phosphatidylcholine{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Choline), C}
const Lysophosphatidylethanolamine{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Ethanolamine), C}
const Phosphatidylethanolamine{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Ethanolamine), C}
const LysophosphatidylNmodifiedethanolamine{B, C} = Monoradylglycerophosphate{B, C} where {D, B <: glycerophospho(Ethanolamine, D), C}
const PhosphatidylNmodifiedethanolamine{B, C} = Diradylglycerophosphate{B, C} where {D, B <: glycerophospho(Ethanolamine, D), C}
const LysophosphatidylNmethylethanolamine{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Nmethylethanolamine), C}
const PhosphatidylNmethylethanolamine{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Nmethylethanolamine), C}
const LysophosphatidylNNdimethylethanolamine{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(NNdimethylethanolamine), C}
const PhosphatidylNNdimethylethanolamine{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(NNdimethylethanolamine), C}
const Lysophosphatidylserine{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Serine), C}
const Phosphatidylserine{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Serine), C}
const LysophosphatidylNmodifiedserine{B, C} = Monoradylglycerophosphate{B, C} where {D, B <: glycerophospho(Serine, D), C}
const PhosphatidylNmodifiedserine{B, C} = Diradylglycerophosphate{B, C} where {D, B <: glycerophospho(Serine, D), C}
const Lysophosphatidylinositol{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Inositol), C}
const Phosphatidylinositol{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Inositol), C} # include PIP
const Lysophosphatidylglycerol{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Glycerol), C}
const Phosphatidylglycerol{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Glycerol), C}
const Lysophosphatidylglycerolphosphate{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophosphophosphate(Glycerol), C}
const Phosphatidylglycerolphosphate{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophosphophosphate(Glycerol), C}
const Lysophosphatidylmethanol{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Methanol), C}
const Phosphatidylmethanol{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Methanol), C}
const Lysophosphatidylethanol{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Ethanol), C}
const Phosphatidylethanol{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Ethanol), C}
const Lysophosphatidylglycan{B, C} = Monoradylglycerophosphate{B, C} where {B <: glycerophospho(Glycan), C}
const Phosphatidylglycan{B, C} = Diradylglycerophosphate{B, C} where {B <: glycerophospho(Glycan), C} # Glcan-GP

struct Bisradylglycerophosphate{B <: DehydratedChemical{<: Tuple{<: includeSIL(Glycerol), <: PPA, <: includeSIL(Glycerol)}}, C} <: Glycerophospholipid{B, C}
    backbone::B
    chain::C
    sn::UInt16
end
const Bisphosphatidicacid{B, C} = Bisradylglycerophosphate{B, C} where {B, C <: TetraradylChain}
# BPA [0] [0, 0] [0, 0, 0] [0, 0, 0, 0] [1/2, 0] [1, PXP, 3] [1, PXP, 3, 3] [1, PXP, 4, 3] [1, 1, PXP, 3, 3] [1, 1, PXP, 4, 3] [1, 2, PXP, 4, 3]
const Semilysobisphosphatidicacid{B, C} = Bisradylglycerophosphate{B, C} where {B, C <: TriradylChain}
# SLBPA [0] [0, 0, 0] [1, 0] [1, P, 3/4] [1/2, 0, P, 0] [1, 1, P, 0] [1, 1, P, 3/4] [1, 2, P, 0] [1, 2, P, 3/4]  
const Lysobisphosphatidicacid{B, C} = Bisradylglycerophosphate{B, C} where {B, C <: DiradylChain} 
# BMP/LBPA [0] [0, 0] [1/2, P, 0] [0, P, 1/2] [1/2, P, 1/2]

struct Bisradylglycerophosphoglycerol{B <: DehydratedChemical{<: Tuple{<: includeSIL(Glycerol), <: PPA, <: includeSIL(Glycerol), <: PPA, <: includeSIL(Glycerol)}}, C} <: Glycerophospholipid{B, C}
    backbone::B
    chain::C
    sn::UInt16
end 
const Cardiolipin{B, C} = Bisradylglycerophosphoglycerol{B, C} where {B, C <: TetraradylChain}
# CL [0] [0, 0] [0, 0, 0] [0, 0, 0, 0] [1/2, 0] [1, PXP, 3] [1, PXP, 3, 3] [1, PXP, 4, 3] [1, 1, PXP, 3, 3] [1, 1, PXP, 4, 3] [1, 2, PXP, 4, 3]
const Monolysocardiolipin{B, C} = Bisradylglycerophosphoglycerol{B, C} where {B, C <: TriradylChain} 
# MLCL [0] [0, 0, 0] [1, 0] [1, PXP, 3/4] [1/2, 0, PXP, 0] [1, 1, PXP, 0] [1, 1, PXP, 3/4] [1, 2, PXP, 0] [1, 2, PXP, 3/4]
const Dilysocardiolipin{B, C} = Bisradylglycerophosphoglycerol{B, C} where {B, C <: DiradylChain} 
# DLCL [0] [0, 0] [1/2, PXP, 0] [0, PXP, 1/2] [1/2, PXP, 1/2]

struct GlycerophosphoNacylethanolamine{B <: glycerophospho(Ethanolamine), C <: CarbonChain{Acyl}} <: Glycerophospholipid{B, C}
    backbone::B
    chain::C
end # GP-NAE
# PnC, PnE, PPA
abstract type Sphingolipid{B, C} <: Lipid{B, C} end
const SUM_CER = Union{<: CarbonChain{Tuple{SPB, Acyl}}, <: Tuple{<: CarbonChain{SPB}, <: CarbonChain{Acyl}}}
const SUM_ACYLCER = Union{<: CarbonChain{Tuple{SPB, Acyl, Acyl}}, <: Tuple{<: CarbonChain{SPB}, <: CarbonChain{Tuple{Acyl, Acyl}}}, <: Tuple{<: CarbonChain{Tuple{SPB, Acyl}}, <: CarbonChain{Acyl}}, <: Tuple{<: CarbonChain{SPB}, <: CarbonChain{Acyl}, <: CarbonChain{Acyl}}}
struct SphingoBone{H, C <: Union{<: CarbonChain{SPB}, <: SUM_CER, <: SUM_ACYLCER}} <: Sphingolipid{H, C}
    headgroup::H
    chain::C
    position::UInt8
end
# HexCer/IPC/SM... headgroup position
# CerP/SPBP.. headgroup position < 32, divrem(x, 32) 35 => [1, 3]
# ACer [0] [0(SPB), 0(ACYL-ACYL)] [0(SPB-ACYL), 0(ACYL)] [2(SPB-ACYL), 0(ACYL)] [1~(SPB-ACYL), 2(ACYL)] [0(SPB), 2(ACYL), 1~(ACYL)] 
# divrem(x, 3) = 
# (0, 2) = [2(SPB-ACYL), 0(ACYL)]
# (a, 2) = [a(SPB-ACYL), 2(ACYL)]/[0(SPB), 2(ACYL), a(ACYL)] 

struct MixSphingoBone{H, C <: Union{<: CarbonChain{SPB}, <: SUM_CER, <: SUM_ACYLCER}} <: Sphingolipid{H, C}
    headgroup::H
    chain::C
    position::Vector{UInt8}
end
# > 1 head group
# pos = (hg1, ..., hg2, chain (optional))
# Hex-(FA....-)ACer/FA....-ASM 
# P-(Hex-)Cer(4, 1) 18:0;3OH/18:0 = Hex-(P-)Cer(1, 4) 18:0;3OH/18:0
# P-(P-)(Hex-)Cer(3, 4, 1) 18:0/18:0 = Hex-(P-)(P-)Cer(1, 3, 4) 18:0/18:0 = P-(Hex-)(P-)Cer(3, 1, 4) 18:0/18:0 = 
# Hex-(P-)(FA 18:0-)ACer(1, 4, 3) 18:0/18:0

const CeramideBone{H, C} = SphingoBone{H, C} where {H, C <: SUM_CER}
const SphingoidBaseBone{H, C} = SphingoBone{H, C} where {H, C <: CarbonChain{SPB}}
const Ceramide{C} = CeramideBone{Nothing, C} where C
const SphingoidBase{C} = SphingoidBaseBone{Nothing, C} where C
const Glycosylceramide{H, C} = CeramideBone{H, C} where {H <: AbstractGlycan, C}
const Glycosylsphingoidbase{H, C} = SphingoidBaseBone{H, C} where {H <: AbstractGlycan, C}

const Manlike = includeSIL(Glycan{<: Tuple{<: Man}})
const Hexlike = includeSIL(Glycan{<: Tuple{<: AbstractHex}})

const CeramidePhosphate{H, C} = CeramideBone{H, C} where {H <: PPA, C}
const SphingoidBasePhosphate{H, C} = SphingoidBaseBone{H, C} where {H <: PPA, C}
const Inositolphosphorylceramide{H, C} = CeramideBone{H, C} where {H <: DehydratedChemical{Tuple{<: PPA, <: includeSIL(Inositol)}}, C}
const Lysoinositolphosphorylceramide{H, C} = SphingoidBaseBone{H, C} where {H <: DehydratedChemical{Tuple{<: PPA, <: includeSIL(Inositol)}}, C}
const Ethanolaminephosphorylceramide{H, C} = CeramideBone{H, C} where {H <: DehydratedChemical{Tuple{<: PPA, <: includeSIL(Ethanolamine)}}, C}
const Lysoethanolaminephosphorylceramide{H, C} = SphingoidBaseBone{H, C} where {H <: DehydratedChemical{Tuple{<: PPA, <: includeSIL(Ethanolamine)}}, C}
const Mannosylinositolphosphorylceramide{H, C} = CeramideBone{H, C} where {H <: DehydratedChemical{<: Tuple{<: PPA, <: includeSIL(Inositol), <: Manlike}}, C}
const Lysomannosylinositolphosphorylceramide{H, C} = SphingoidBaseBone{H, C} where {H <: DehydratedChemical{Tuple{<: PPA, <: includeSIL(Inositol), <: Manlike}}, C}
const Mannosyldiinositolphosphorylceramide{H, C} = CeramideBone{H, C} where {H <: DehydratedChemical{<: Tuple{<: PPA, <: includeSIL(Inositol), <: Manlike, <: PPA, <: includeSIL(Inositol)}}, C}
const Lysomannosyldiinositolphosphorylceramide{H, C} = SphingoidBaseBone{H, C} where {H <: DehydratedChemical{Tuple{<: PPA, <: includeSIL(Inositol), <: Manlike, <: PPA, <: includeSIL(Inositol)}}, C}
const Glycosylinositolphosphorylceramide{H, C} = CeramideBone{H, C} where {H <: DehydratedChemical{Tuple{<: PPA, <: includeSIL(Inositol), Glycan}}, C}
const Lysoglycosylinositolphosphorylceramide{H, C} = SphingoidBaseBone{H, C} where {H <: DehydratedChemical{Tuple{<: PPA, <: includeSIL(Inositol), <: Glycan}}, C}
const Sphingomyelin{H, C} = CeramideBone{H, C} where {H <: DehydratedChemical{Tuple{<: PPA, <: includeSIL(Choline)}}, C}
const Lysosphingomyelin{H, C} = SphingoidBaseBone{H, C} where {H <: DehydratedChemical{Tuple{<: PPA, <: includeSIL(Choline)}}, C}

const Sulfonolipid{H, C} = CeramideBone{H, C} where {H <: includeSIL(SulfurousAcid), C}
const Lysosulfonolipid{H, C} = SphingoidBaseBone{H, C} where {H <: includeSIL(SulfurousAcid), C}

const Acylceramide{C} = SphingoBone{Nothing, C} where {C <: SUM_ACYLCER}
const Acylhexosylceramide{H, C} = MixSphingoBone{H, C} where {H <: Hexlike, C <: SUM_ACYLCER}
const Acylsphingomyelin{H, C} = MixSphingoBone{H, C} where {H <: DehydratedChemical{Tuple{<: PPA, <: includeSIL(Choline)}}, C <: SUM_ACYLCER}

abstract type Sterol{C} <: Lipid{Nothing, C} end
struct SterolBone{C} <: Sterol{C}
    chain::C
end
const FreeSterol{C} = SterolBone{C} where {C <: CarbonChain{<: AbstractSTRing}}
const Sterylester{C} = SterolBone{C} where {C <: Tuple{<: CarbonChain{<: AbstractSTRing}, <: CarbonChain{Acyl}}}

struct SubstitutedSterol{C, S, T} <: Sterol{C}
    backbone::SubstitutedChemical{FreeSterol{C}, S, T}
end
# ST;..., BA;..., SG, ASG

abstract type Prenol{C} <: Lipid{Nothing, C} end
struct Retinylester{C <: CarbonChain{<: Acyl}} <: Prenol{C}
    chain::C
end
struct CoenzymeQ{C <: IsoprenoidChain} <: Prenol{C} end

include("interface.jl")
include("io.jl")
 
end