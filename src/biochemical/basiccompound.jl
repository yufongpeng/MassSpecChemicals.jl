module BasicCompounds
using Reexport
@reexport using ..BioChemicals
using ..MSChemicals: AbstractChemical
export BasicCompound, 
       Dihydrogen,
       Hydrogen,
       Ethane, 
       Ethyl, # Et
       Ethanol, 
       Ethoxy, # OEt
       Methane, 
       Methyl, # Me
       Methanol, 
       Methoxy, # OMe
       FormicAcid, 
       Oformyl, # OFo
       Formamide, 
       Nformyl, # NFo
       AceticAcid, 
       Oacetyl, # OAc
       Acetamide, 
       Nacetyl, # NAc
       HydrogenCyanide,
       Cyano, # CN

       HydrogenBromide, 
       Bromo, # Br
       HydrogenChloride, 
       Chloro, # Cl
       HydrogenFluoride, 
       Fluoro, # F
       HydrogenIodide, 
       Iodo, # I
       HydrogenOxide, 
       OxygenAtom, 
       Hydroxy, # OH
       CarboxylicAcidGroup, # COOH
       Oxo, # oxo 2
       Alkoxy, # oxy 2
       HydrogenPeroxide, 
       Hydroperoxyl, # OOH
       Epoxy, # Ep 2
       Peroxy, # OO 2
       Ammonia, 
       Amino, # NH2/N
       HydrogenSulfide, 
       Sulfanyl, # SH
       HydrogenCyanide, 
       Cyano, # CN
       PhosphoricAcid, # P
       Phosphate, # P/MP
       DiphosphoricAcid, 
       Diphosphate, # DP
       TriphosphoricAcid, 
       Triphosphate, # TP
       SulfuricAcid, # S
       Sulfate, # S
       SulfurousAcid, 
       Sulfite,
       Formyl, # Fo
       Acetyl, # Ac
       NitricAcid, 
       Nitro, # NO2
       Glycerol

abstract type BasicCompound <: AbstractChemical end 

struct Dihydrogen <: BasicCompound end 
struct Hydrogen <: FunctionalGroup{Dihydrogen, Dehydrogen} end # H
struct Ethane <: BasicCompound end 
struct Ethyl <: FunctionalGroup{Ethane, Dehydrogen} end # Et
struct Ethanol <: BasicCompound end 
struct Ethoxy <: FunctionalGroup{Ethanol, Dehydrogen} end # OEt
struct Methane <: BasicCompound end 
struct Methyl <: FunctionalGroup{Methane, Dehydrogen} end # Me
struct Methanol <: BasicCompound end 
struct Methoxy <: FunctionalGroup{Methanol, Dehydrogen} end # OMe
struct FormicAcid <: BasicCompound end 
struct Oformyl <: FunctionalGroup{FormicAcid, Dehydrogen} end # OFo
struct Formamide <: BasicCompound end 
struct Nformyl <: FunctionalGroup{Formamide, Dehydrogen} end # NFo
struct AceticAcid <: BasicCompound end 
struct Oacetyl <: FunctionalGroup{AceticAcid, Dehydrogen} end # OAc
struct Acetamide <: BasicCompound end 
struct Nacetyl <: FunctionalGroup{Acetamide, Dehydrogen} end # NAc
struct HydrogenCyanide <: BasicCompound end
struct Cyano <: FunctionalGroup{HydrogenCyanide, Dehydrogen} end # CN
struct HydrogenBromide <: BasicCompound end
struct Bromo <: FunctionalGroup{HydrogenBromide, Dehydrogen} end # Br
struct HydrogenChloride <: BasicCompound end
struct Chloro <: FunctionalGroup{HydrogenChloride, Dehydrogen} end # Cl
struct HydrogenFluoride <: BasicCompound end
struct Fluoro <: FunctionalGroup{HydrogenFluoride, Dehydrogen} end # F
struct HydrogenIodide <: BasicCompound end
struct Iodo <: FunctionalGroup{HydrogenIodide, Dehydrogen} end # I
struct HydrogenOxide <: BasicCompound end
struct OxygenAtom <: FunctionalGroup{HydrogenOxide, Didehydrogen} end 
struct Hydroxy <: FunctionalGroup{HydrogenOxide, Dehydrogen} end # OH
struct CarboxylicAcidGroup <: FunctionalGroup{FormicAcid, Dehydrogen} end # COOH
struct Oxo <: FunctionalGroup{HydrogenOxide, Didehydrogen} end # oxo 2
struct Alkoxy <: FunctionalGroup{HydrogenOxide, Didehydrogen} end # oxy 2
struct HydrogenPeroxide <: BasicCompound end
struct Hydroperoxyl <: FunctionalGroup{HydrogenPeroxide, Dehydrogen} end # OOH
struct Epoxy <: FunctionalGroup{HydrogenOxide, Didehydrogen} end # Ep 2
struct Peroxy <: FunctionalGroup{HydrogenPeroxide, Didehydrogen} end # OO 2
struct Ammonia <: BasicCompound end
struct Amino <: FunctionalGroup{Ammonia, Dehydrogen} end # NH2/N
struct HydrogenSulfide <: BasicCompound end
struct Sulfanyl <: FunctionalGroup{HydrogenSulfide, Dehydrogen} end # SH
struct PhosphoricAcid <: BasicCompound end # P
struct Phosphate <: FunctionalGroup{PhosphoricAcid, Dehydrogen} end # P/MP
struct DiphosphoricAcid <: BasicCompound end
struct Diphosphate <: FunctionalGroup{DiphosphoricAcid, Dehydrogen} end # DP
struct TriphosphoricAcid <: BasicCompound end
struct Triphosphate <: FunctionalGroup{TriphosphoricAcid, Dehydrogen} end # TP
struct SulfuricAcid <: BasicCompound end # S
struct Sulfate <: FunctionalGroup{SulfuricAcid, Dehydrogen} end # S
struct SulfurousAcid <: BasicCompound end
struct Sulfite <: FunctionalGroup{SulfurousAcid, Dehydrogen} end
struct Formyl <: FunctionalGroup{FormicAcid, Dehydroxy} end # Fo
struct Acetyl <: FunctionalGroup{AceticAcid, Dehydroxy} end # Ac
struct NitricAcid <: BasicCompound end
struct Nitro <: FunctionalGroup{NitricAcid, Dehydroxy} end # NO2
struct Glycerol <: BasicCompound end
end