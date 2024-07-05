const AA_3LETTER = ["Gly", "Ser", "Tyr", "Tau"]
const GLYCOLIPID = ["GM3", "GM2", "GM1a", "GD1a", "GD1aα", "GD1aa", "GT1a", "GT1aα", "GT1aa", "GD3", "GD2", "GD1b", "GT1b", "GT1bα", "GT1ba", "GQ1b", "GQ1bα", "GQ1ba",
"GT3", "GT2", "GT1c", "GQ1c", "GQ1cα", "GQ1ca", "GP1c", "GP1cα", "GP1ca", "GA2", "GA1", "GM1b", "GM1α", "GM1alpha", "GD1c", "GD1α", "GD1alpha", "GD1e",
"Gb3", "Gb4", "iGb3", "iGb4", "Lc3", "LM1", "GM4", "SM4", "SM3", "SM2", "SM1a", "SM1b", "SB1a"]
const FREESTEROL = ["C", "DSMS", "DC", "CAS", "BRS", "EGS", "DEGS", "SIS", "STS"]
const CLASS = [
    "HC", "FA", "FAL", "FOH", "WE", "FAM", "FattyAmide", "NA", "NAE", [string("NA", x) for x in AA_3LETTER]..., "CAR", "CoA", "FAHFA",
    "MG", "DG", "TG", "SQMG", "SQDG", "MGMG", "MGDG", "DGMG", "DGDG", "GlcAMG", "GlcADG",
    "PA", "LPA", "PC", "LPC", "PE", "LPE", "PE-NMe", "LPE-NMe", "PE-NMe2", "LPE-NMe2", "PE-N", "LPE-N", "PS", "LPS", "PS-N", "LPS-N", "PI", "LPI",
    "PG", "LPG", "PMeOH", "LPMeOH", "PEtOH", "LPEtOH", "GP", "LGP", "PGP", "LPGP", "BPA", "SLBPA", "BMP", "LBPA", "CL", "MLCL", "DLCL", "GP-NAE", "PIP", "PIP2", "PIP3",
    "Cer", "SPB", "CerP", "SPBP", "EPC", "LEPC", "PE-Cer", "PE-SPB", "IPC", "LIPC", "PI-Cer", "PI-SPB", "MIPC", "LMIPC", "M(IP)2C", "LM(IP)2C", "SM", "LSM",
    GLYCOLIPID..., [string("Lyso", x) for x in GLYCOLIPID]..., "SL", "LSL", "ACer", "ASM", 
    "ST", "SE", "BA", "SG", "ASG", (FREESTEROL .* "E")..., "FC", FREESTEROL[begin + 1:end]...
]

const POST_MODIFIER = Dict{String, String}(
    "LPE-N" => "(\\([^)(]*+(?:(?1)[^)(]*)*+\\))",
    "LPS-N" => "(\\([^)(]*+(?:(?1)[^)(]*)*+\\))",
    "PE-N"  => "(\\([^)(]*+(?:(?1)[^)(]*)*+\\))",
    "PS-N"  => "(\\([^)(]*+(?:(?1)[^)(]*)*+\\))"
)

const PRE_MODIFIER = Dict{String, String}(
    "ASM"       => "FA\\s+\\d+:\\d+\\s*(?:\\([\\d,\\,,E,Z,\\s]+\\))?\\s*(?:\\[[^\\]\\[]+\\])?((?:;[^;\\(/_]*(?:\\([^)(]*+(?:(?1)[^)(]*)*+\\))?[^;/_]*)*)\\-",
    "ACer"      => "FA\\s+\\d+:\\d+\\s*(?:\\([\\d,\\,,E,Z,\\s]+\\))?\\s*(?:\\[[^\\]\\[]+\\])?((?:;[^;\\(/_]*(?:\\([^)(]*+(?:(?1)[^)(]*)*+\\))?[^;/_]*)*)\\-",
    "AHexCer"   => "FA\\s+\\d+:\\d+\\s*(?:\\([\\d,\\,,E,Z,\\s]+\\))?\\s*(?:\\[[^\\]\\[]+\\])?((?:;[^;\\(/_]*(?:\\([^)(]*+(?:(?1)[^)(]*)*+\\))?[^;/_]*)*)\\-",
    "AGlcCer"   => "FA\\s+\\d+:\\d+\\s*(?:\\([\\d,\\,,E,Z,\\s]+\\))?\\s*(?:\\[[^\\]\\[]+\\])?((?:;[^;\\(/_]*(?:\\([^)(]*+(?:(?1)[^)(]*)*+\\))?[^;/_]*)*)\\-",
    "AGalCer"   => "FA\\s+\\d+:\\d+\\s*(?:\\([\\d,\\,,E,Z,\\s]+\\))?\\s*(?:\\[[^\\]\\[]+\\])?((?:;[^;\\(/_]*(?:\\([^)(]*+(?:(?1)[^)(]*)*+\\))?[^;/_]*)*)\\-",
)

const CLASS_STRUCT = Dict{String, Function}(
    "HC"    => (mod, cls, pos, sil) -> (FattyAcyl, makemolecule(Dihydrogen(); sil), Alkyl), 
    "FA"    => (mod, cls, pos, sil) -> (FattyAcyl, makemolecule(HydrogenOxide(); sil), Acyl), 
    "FAL"   => (mod, cls, pos, sil) -> (FattyAcyl, makemolecule(Dihydrogen(); sil), Acyl), 
    "FOH"   => (mod, cls, pos, sil) -> (FattyAcyl, makemolecule(HydrogenOxide(); sil), Alkyl), 
    "WE"    => (mod, cls, pos, sil) -> (FattyAcyl, makemolecule(HydrogenOxide(); sil), (Alkyl, Acyl)), 
    "FAM"   => (mod, cls, pos, sil) -> (FattyAcyl, makemolecule(Ammonia(); sil), Acyl), 
    "FattyAmide" => (mod, cls, pos, sil) -> (FattyAcyl, makemolecule(Ammonia(); sil), Acyl), 
    "NA"    => (mod, cls, pos, sil) -> (FattyAcyl, makemolecule(Ammonia(); sil), (Alkyl, Acyl)), 
    "NAE"   => (mod, cls, pos, sil) -> (FattyAcyl, makemolecule(Ethanolamine(); sil), Acyl),
    [string("NA", x) => (mod, cls, pos, sil) -> (FattyAcyl, makemolecule(parse_aa(replace(cls, "NA" => "")); sil), Acyl) for x in AA_3LETTER]...,
    "CAR"   => (mod, cls, pos, sil) -> (FattyAcyl, makemolecule(Carnitine(); sil), Acyl),
    "CoA"   => (mod, cls, pos, sil) -> (FattyAcyl, makemolecule(CoA(); sil), Acyl),
    "FAHFA" => (mod, cls, pos, sil) -> (FattyAcyl, makemolecule(HydrogenOxide(); sil), (Acyl, Acyl)),
    "MG"    => (mod, cls, pos, sil) -> (isnothing(mod) ? Monoradylglycerol : Omodifiedradylglycerol, makemolecule(DehydratedChemical, makemolecule(parse_headgroup(mod)), Glycerol(); sil), Radyl),
    "DG"    => (mod, cls, pos, sil) -> (isnothing(mod) ? Diradylglycerol : Omodifiedradylglycerol, makemolecule(DehydratedChemical, makemolecule(parse_headgroup(mod)), Glycerol(); sil), (Radyl, Radyl)),
    "TG"    => (mod, cls, pos, sil) -> (Triradylglycerol, makemolecule(Glycerol(); sil), (Radyl, Radyl, Radyl)),
    "SQMG"  => (mod, cls, pos, sil) -> (Omodifiedradylglycerol, makemolecule(DehydratedChemical, Sulfoquinovose(), Glycerol(); sil, linkage = [α(0x01) => lk(0x03)]), Radyl),
    "SQDG"  => (mod, cls, pos, sil) -> (Omodifiedradylglycerol, makemolecule(DehydratedChemical, Sulfoquinovose(), Glycerol(); sil, linkage = [α(0x01) => lk(0x03)]), (Radyl, Radyl)),
    "MGMG"  => (mod, cls, pos, sil) -> (Omodifiedradylglycerol, makemolecule(DehydratedChemical, Gal(), Glycerol(); sil, linkage =  [β(0x01) => lk(0x03)]), Radyl),
    "MGDG"  => (mod, cls, pos, sil) -> (Omodifiedradylglycerol, makemolecule(DehydratedChemical, Gal(), Glycerol(); sil, linkage =  [β(0x01) => lk(0x03)]), (Radyl, Radyl)),
    "DGMG"  => (mod, cls, pos, sil) -> (Omodifiedradylglycerol, makemolecule(DehydratedChemical, Gal(), Gal(), Glycerol(); sil, linkage =  [α(0x01) => lk(0x06), β(0x01) => lk(0x03)]), Radyl),
    "DGDG"  => (mod, cls, pos, sil) -> (Omodifiedradylglycerol, makemolecule(DehydratedChemical, Gal(), Gal(), Glycerol(); sil, linkage =  [α(0x01) => lk(0x06), β(0x01) => lk(0x03)]), (Radyl, Radyl)),
    "GlcAMG"    => (mod, cls, pos, sil) -> (Omodifiedradylglycerol, makemolecule(DehydratedChemical, GlcA(), Glycerol(); sil, linkage =  [α(0x01) => lk(0x03)]), Radyl),
    "GlcADG"    => (mod, cls, pos, sil) -> (Omodifiedradylglycerol, makemolecule(DehydratedChemical, GlcA(), Glycerol(); sil, linkage =  [α(0x01) => lk(0x03)]), (Radyl, Radyl)),
    "PA"    => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "LPA"   => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, PhosphoricAcid(), Glycerol(); sil), Radyl, Radyl), 
    "PC"    => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Choline(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "LPC"   => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Choline(), PhosphoricAcid(), Glycerol(); sil), Radyl), 
    "PE"    => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Ethanolamine(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "LPE"   => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Ethanolamine(), PhosphoricAcid(), Glycerol(); sil), Radyl), 
    "PE-NMe"    => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Nmethylethanolamine(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "LPE-NMe"   => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Nmethylethanolamine(), PhosphoricAcid(), Glycerol(); sil), Radyl), 
    "PE-NMe2"   => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, NNdimethylethanolamine(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "LPE-NMe2"  => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, NNdimethylethanolamine(), PhosphoricAcid(), Glycerol(); sil), Radyl), 
    "PE-N"  => (mod, cls, pos, sil) -> (Radylglycerophosphate, concatmolecule(DehydratedChemical, makemolecule(parse_Nmodification(cls)), makemolecule(Ethylethanolamine(), PhosphoricAcid(), Glycerol(); sil)), (Radyl, Radyl)), 
    "LPE-N" => (mod, cls, pos, sil) -> (Radylglycerophosphate, concatmolecule(DehydratedChemical, makemolecule(parse_Nmodification(cls)), makemolecule(Ethylethanolamine(), PhosphoricAcid(), Glycerol(); sil)), Radyl), 
    "PS"    => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Serine(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "LPS"   => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Serine(), PhosphoricAcid(), Glycerol(); sil), Radyl), 
    "PS-N"  => (mod, cls, pos, sil) -> (Radylglycerophosphate, concatmolecule(DehydratedChemical, makemolecule(parse_Nmodification(cls)), makemolecule(Serine(), PhosphoricAcid(), Glycerol(); sil)), (Radyl, Radyl)), 
    "LPS-N" => (mod, cls, pos, sil) -> (Radylglycerophosphate, concatmolecule(DehydratedChemical, makemolecule(parse_Nmodification(cls)), makemolecule(Serine(), PhosphoricAcid(), Glycerol(); sil)), Radyl), 
    "PI"    => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Inositol(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "LPI"   => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Inositol(), PhosphoricAcid(), Glycerol(); sil), Radyl), 
    "PG"    => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Glycerol(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "LPG"   => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Glycerol(), PhosphoricAcid(), Glycerol(); sil), Radyl), 
    "PMeOH" => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Methanol(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "LPMeOH"    => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Methanol(), PhosphoricAcid(), Glycerol(); sil), Radyl), 
    "PEtOH" => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Ethanol(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "LPEtOH"    => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, Ethanol(), PhosphoricAcid(), Glycerol(); sil), Radyl), 
    "GP"    => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, makemolecule(parse_headgroup(mod)), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "LGP"   => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, makemolecule(parse_headgroup(mod)), PhosphoricAcid(), Glycerol(); sil), Radyl), 
    "PGP"   => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, PhosphoricAcid(), Glycerol(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "LPGP"  => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, PhosphoricAcid(), Glycerol(), PhosphoricAcid(), Glycerol(); sil), Radyl), 
    "BPA"   => (mod, cls, pos, sil) -> (Bisradylglycerophosphate, makemolecule(DehydratedChemical, Glycerol(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl, Radyl, Radyl)), 
    "SLBPA" => (mod, cls, pos, sil) -> (Bisradylglycerophosphate, makemolecule(DehydratedChemical, Glycerol(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl, Radyl)), 
    "BMP"   => (mod, cls, pos, sil) -> (Bisradylglycerophosphate, makemolecule(DehydratedChemical, Glycerol(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "LBPA"  => (mod, cls, pos, sil) -> (Bisradylglycerophosphate, makemolecule(DehydratedChemical, Glycerol(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "CL"    => (mod, cls, pos, sil) -> (Bisradylglycerophosphoglycerol, makemolecule(DehydratedChemical, Glycerol(), PhosphoricAcid(), Glycerol(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl, Radyl, Radyl)), 
    "MLCL"  => (mod, cls, pos, sil) -> (Bisradylglycerophosphoglycerol, makemolecule(DehydratedChemical, Glycerol(), PhosphoricAcid(), Glycerol(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl, Radyl)), 
    "DLCL"  => (mod, cls, pos, sil) -> (Bisradylglycerophosphoglycerol, makemolecule(DehydratedChemical, Glycerol(), PhosphoricAcid(), Glycerol(), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "GP-NAE"    => (mod, cls, pos, sil) -> (GlycerophosphoNacylethanolamine, makemolecule(DehydratedChemical, Ethanolamine(), PhosphoricAcid(), Glycerol(); sil), Acyl), 
    "PIP"   => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, isnothing(pos) ? Inositol([Phosphate()]) : Inositol([parse(UInt8, x.match) => Phosphate() for x in eachmatch(r"\d+", pos)][begin:begin]), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "PIP2"  => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, isnothing(pos) ? Inositol([Phosphate(), Phosphate()]) : Inositol([parse(UInt8, x.match) => Phosphate() for x in eachmatch(r"\d+", pos)][begin:begin + 1]), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)), 
    "PIP3"  => (mod, cls, pos, sil) -> (Radylglycerophosphate, makemolecule(DehydratedChemical, isnothing(pos) ? Inositol([Phosphate(), Phosphate(), Phosphate()]) : Inositol([parse(UInt8, x.match) => Phosphate() for x in eachmatch(r"\d+", pos)][begin:begin + 2]), PhosphoricAcid(), Glycerol(); sil), (Radyl, Radyl)),
    "Cer"   => (mod, cls, pos, sil) -> isnothing(mod) ? (SphingoBone, nothing, (SPB, Acyl)) : endswith(mod, ")") ? (MixedSphingoBone, makemolecule(parse_headgroup(mod); sil), (SPB, Acyl)) : (SphingoBone, makemolecule(parse_headgroup(mod); sil), (SPB, Acyl)),
    "SPB"   => (mod, cls, pos, sil) -> isnothing(mod) ? (SphingoBone, nothing, SPB) : endswith(mod, ")") ? (MixedSphingoBone, makemolecule(parse_headgroup(mod); sil), SPB) : (SphingoBone, makemolecule(parse_headgroup(mod); sil), SPB),
    "CerP"  => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(PhosphoricAcid(); sil), (SPB, Acyl)),
    "SPBP"  => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(PhosphoricAcid(); sil), SPB),
    "EPC"   => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(DehydratedChemical, Ethanolamine(), PhosphoricAcid(); sil), (SPB, Acyl)),
    "LEPC"  => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(DehydratedChemical, Ethanolamine(), PhosphoricAcid(); sil), SPB),
    "PE-Cer"    => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(DehydratedChemical, Ethanolamine(), PhosphoricAcid(); sil), (SPB, Acyl)),
    "PE-SPB"    => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(DehydratedChemical, Ethanolamine(), PhosphoricAcid(); sil), SPB),
    "IPC"   => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(DehydratedChemical, Inositol(), PhosphoricAcid(); sil), (SPB, Acyl)),
    "LIPC"  => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(DehydratedChemical, Inositol(), PhosphoricAcid(); sil), SPB),
    "PI-Cer"    => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(DehydratedChemical, Inositol(), PhosphoricAcid(); sil), (SPB, Acyl)),
    "PI-SPB"    => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(DehydratedChemical, Inositol(), PhosphoricAcid(); sil), SPB),
    "MIPC"  => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(DehydratedChemical, Man(), Inositol(), PhosphoricAcid(); sil, linkage = [α(0x01) => lk(0x02), α(0x01) => lk(0x00)]), (SPB, Acyl)),
    "LMIPC" => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(DehydratedChemical, Inositol(), Man(), PhosphoricAcid(); sil, linkage = [α(0x01) => lk(0x02), α(0x01) => lk(0x00)]), SPB),
    "M(IP)2C"   => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(DehydratedChemical, Inositol(), PhosphoricAcid(), Man(), Inositol(), PhosphoricAcid(); sil, linkage = [α(0x01) => lk(0x00), lk(0x00) => lk(0x00), α(0x01) => lk(0x02), α(0x01) => lk(0x00)]), (SPB, Acyl)),
    "LM(IP)2C"  => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(DehydratedChemical, Inositol(), PhosphoricAcid(), Man(), Inositol(), PhosphoricAcid(); sil, linkage = [α(0x01) => lk(0x00), lk(0x00) => lk(0x00), α(0x01) => lk(0x02), α(0x01) => lk(0x00)]), SPB),
    "SM"    => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(DehydratedChemical, Choline(), PhosphoricAcid(); sil), (SPB, Acyl)),
    "LSM"   => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(DehydratedChemical, Choline(), PhosphoricAcid(); sil), SPB),
    [cls => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(eval(cls)(); sil), (SPB, Acyl)) for cls in GLYCOLIPID]...,
    [string("Lyso", cls) => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(eval(cls)(); sil), SPB) for cls in GLYCOLIPID]...,
    "SL"    => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(SulfurousAcid(); sil), (SPB, Acyl)),
    "LSL"   => (mod, cls, pos, sil) -> (SphingoBone, makemolecule(SulfurousAcid(); sil), SPB),
    "ACer"  => (mod, cls, pos, sil) -> (SphingoBone, [parse_acylcerchain(Acyl, cls; position = pos)], (SPB, Acyl, Acyl)),
    "ASM"   => (mod, cls, pos, sil) -> (MixSphingoBone, [makemolecule(DehydratedChemical, Choline(), PhosphoricAcid(); sil), parse_acylcerchain(Acyl, cls; position = pos)], (SPB, Acyl, Acyl)),
    # AHexCer ?
    "ST"    => (mod, cls, pos, sil) -> (SterolBone, nothing, STRing),
    "SE"    => (mod, cls, pos, sil) -> (SterolBone, nothing, (STRing, Acyl)),
    "SG"    => (mod, cls, pos, sil) -> (SubstitutedSterol, nothing, STRing),
    "ASG"   => (mod, cls, pos, sil) -> (SubstitutedSterol, nothing, STRing),
    "FC"    => (mod, cls, pos, sil) -> (SterolBone, nothing, CRing),
    [cls => (mod, cls, pos, sil) -> (SterolBone, nothing, eval(string(cls, "Ring"))) for cls in FREESTEROL[begin + 1:end]]..., 
    [string(cls, "E") => (mod, cls, pos, sil) -> (SterolBone, nothing, (eval(string(cls, "Ring")), Acyl)) for cls in FREESTEROL]..., 
    [string(cls, "G") => (mod, cls, pos, sil) -> (SubstitutedSterol, nothing, eval(string(cls, "Ring"))) for cls in FREESTEROL]..., 
    [string("A", cls, "G") => (mod, cls, pos, sil) -> (SubstitutedSterol, nothing, eval(string(cls, "Ring"))) for cls in FREESTEROL]..., 
    "BA"    => (mod, cls, pos, sil) -> (SterolBone, nothing, BARing)
)

function class_regex()
    C = collect(keys(CLASS_STRUCT))
    OC = map(C) do x
        replace(x, "(" => "\\(", ")" => "\\)", "-" => "\\-")
    end
    RC = map(C) do x                                                                                                                                   
        get(PRE_MODIFIER, x, "") * replace(x, "(" => "\\(", ")" => "\\)", "-" => "\\-") * get(POST_MODIFIER, x, "")                                                                                                                                              
    end
    sort!(RC; by = x -> length(x), rev = true)
    # sort!(C; by = x -> length(x), rev = true)
    KC = join(RC, "|")
    DC = join(OC, ")(?:")
    Regex("^(.*?)(" * KC * ")[^(?:" * DC * ")]*(\\([\\d,\\,,',\\s]+\\))?\\s*(\\[[^\\]\\[]+\\])?[\\s,;]")
end

function chain_regex()
    r"([\s,/,_][d,t,e]?[P,O]?-?\d+:\d+)(\([^)(]*+(?:(?1)[^)(]*)*+\))?(\[[^\[]+\])?((?:;[^)(;(/_]*(\([^)(]*+(?:(?2)[^)(]*)*+\))?[^)(;(/_]*)*)(\(sn-*\d*'*\))?"
end

function clsmod_regex()
    pre = collect(keys(PRE_MODIFIER))
    post = collect(keys(POST_MODIFIER))
    C = [pre..., post...]
    unique!(C)
    RC = map(C) do x
        replace(x, "(" => "\\(", ")" => "\\)", "-" => "\\-")
    end
    sort!(RC; by = x -> length(x), rev = true)
    KC = join(RC, "|")
    DC = join(RC, ")(?:")
    Regex("^(.*?)(" * KC * ")[^(?:" * DC * ")]*")
end

const REGEX = Dict{Symbol, Regex}(
    :class => class_regex(),
    :chain => chain_regex(),
    :clsmod => clsmod_regex()
)

function register_class!(cls::String)
    push!(CLASS, cls)
    REGEX[:class] = class_regex()
end

const SPECIFIC_LIPID =  ["Cholesterol", "Desmosterol", "Dihydrocholesterol", "Campesterol", "Brassicasterol", "Ergosterol", "Dehydroergosterol", "Sitosterol", "Stigmasterol"]
function parse_lipid(s::AbstractString)
    any(==(s), SPECIFIC_LIPID) && return parse_spesific(s)
    result = match(REGEX[:class], s)
    mod = result.captures[begin]
    mod = isnothing(mod) ? nothing : isempty(mod) ? nothing : mod
    cls = result.captures[begin + 1]
    pos = result.captures[end - 1]
    pos = isnothing(pos) ? nothing : isempty(pos) ? nothing : pos
    sil = result.captures[end]
    sil = isnothing(sil) ? nothing : isempty(sil) ? nothing : sil
    parse_head = get(CLASS_STRUCT, string(cls), nothing)
    parse_head = isnothing(parse_head) ? get(CLASS_STRUCT, string(match(REGEX[:clsmod], cls).captures[2]), nothing) : parse_head
    # parse head first, parse chain with head
    Con, Bone, Chain = parse_head(mod, cls, pos, sil)
    schain = s[result.match.ncodeunits:end]
    parse_chain(Con, Bone, Chain, schain)
end
# match(r"([\s,/,_][d,t,e]?[P,O]?-?\d+:\d+)(\([^)(]*+(?:(?1)[^)(]*)*+\))?(\[[^\[]+\])?((?:;[^;\(/_]*(?:\([^)(]*+(?:(?1)[^)(]*)*+\))?[^;/_]*)*)", s)
# cbdb, pos, fg = match(r"[/,_]?([d,t,e]?[P,O]?-?\d*:\d*)(\([^)(]*+(?:(?1)[^)(]*)*+\))?(\[[^\[]+\])?((?:;[^;\(/_]*(?:\([^)(]*+(?:(?1)[^)(]*)*+\))?[^;/_]*)*)", s)
function parse_chain(Con::Type{FattyAcyl}, Bone, Chain, schain)
    mchain = collect(eachmatch(REGEX[:chain], schain))
    if length(Chain) == mchain
        [parse_fattychain(C, m) for (C, m) in zip(Chain, mchain)]
    else
        parse_fattychain(Chain, mchain)
    end
end

function parse_chain(Con::Type{<: Union{<: Glycerolipid, <: Glycerophospholipid}}, Bone, Chain, schain)
    mchain = collect(eachmatch(REGEX[:chain], schain))
    maxsn = length(snposition(Con))
    pchain = if length(mchain) == maxsn
        [parse_fattychain(Radyl, m) for m in mchain] # ordered/unordered
    elseif length(mchain) == length(Chain) # unordered
        [parse_fattychain(Radyl, m) for m in mchain]
    else
        snx = length(Chain) ÷ length(mchain)
        sns = repeat([snx], length(mchain))
        δ = length(Chain) - length(mchain) * snx
        i = lastindex(sns)
        while δ > 0
            δ -= 1
            sns[i] += 1
            i -= 1
        end
        [sn == 1 ? parse_fattychain(Radyl, m) : parse_fattychain(ntuple(i -> Radyl, sn), m) for (sn, m) in zip(sns, mchain)]
    end
    # check sn position distribution (sn)
    # :s => 0x01
    # :o => vector ID 
    # :u => 0x00
    # sn => snposition ID 
    # base = maxsn + 1
    position = map(enumerate(pchain)) do (i, c)
        @match first(c) begin 
            :s => 0x01
            :o => UInt8(i)
            :u => 0x00
            r"\(sn.*\)" => UInt8(findfirst(==(replace(c, "(" => "", ")" => "")), snposition(Con)))
        end
    end
    sn = maxsn > 3 ? 0x0000 : 0x00
    bs = convert(typeof(sn), maxsn + 1)
    for p in position
        sn *= bs
        sn += p
    end
    Con(Bone, length(pchain) == 1 ? last(first(pchain)) : ntuple(i -> last(pchain[i]), length(pchain)), sn)
    # prev = prev * base + next
    # check sn chain number == maxsn for sn-position level
end
function parse_chain(Con::Type{<: Sphingolipid}, Bone, Chain, schain)
    mchain = collect(eachmatch(REGEX[:chain], schain))
    maxsn = length(snposition(Con))
    if length(Chain) == 3 # ACer...
        # take chain from Bone
        addchain = pop!(Bone)
        if isnothing(addchain) && length(mchain) == 1 # sum level
            pchain = parse_fattychain(Chain, mchain)
        elseif length(mchain) == 1 # known additional acyl on lcb
            pchain = [parse_fattychain(Chain[begin:begin + 1], mchain)]
            push!(pchain, addchain)
        elseif length(mchain) == 2 # species level
            pchain = [parse_fattychain(C, m) for (C, m) in zip(Chain[begin:begin + 1], mchain)]
            push!(pchain, addchain)
        end
    elseif length(mchain) == length(Chain)
        pchain = [parse_fattychain(C, m) for (C, m) in zip(Chain, mchain)]
    else
        pchain = parse_fattychain(Chain, mchain)
    end
    position = map(enumerate(pchain)) do (i, c)
        @match first(c) begin 
            :s => 0x01
            :o => UInt8(i)
        end
    end
    sn = 0x00
    bs = convert(typeof(sn), maxsn + 1)
    for p in position
        sn *= bs
        sn += p
    end
    Con(Bone, length(pchain) == 1 ? last(first(pchain)) : ntuple(i -> last(pchain[i]), length(pchain)), sn)
end

function parse_fattychain(T, mchain::RegexMatch; position = nothing)
    # pos => chain
    # pos:
    # :s start
    # :o ordered
    # :u unordered
    # number for ACer
    # sn-\d
    cchain, pos, sil, mod, _, sn = mchain
    p = !isnothing(position) ? position : 
        !isnothing(sn) ? sn : 
        startswith(mchain.match, r"\s") ? :s : 
        startswith(mchain.match, "/") ? :o :
        startswith(mchain.match, "_") ? :u : throw(ArgumentError("Invalid fattyacyl chain."))
    p => parse_carbonchain(T, cchain, pos, sil, mod)
end

function parse_carbonchain(T::Tuple, cchain, pos, sil, mod)
    if SPB in T
        CarbonChain{Tuple{T...}}(parse_singlechain(SPB, cchain, pos, sil, mod)...)
    elseif any(x -> supertype(x) == AbstractSTRing, T)
        i = findfirst(x -> supertype(x) == AbstractSTRing, T)
        CarbonChain{Tuple{T...}}(parse_singlechain(T[i], cchain, pos, sil, mod)...)
    else
        rad, cchain = match(r"([d,t,e]?[P,O]?-?)(\d+:\d+.*)", cchain)
        if isnothing(rad) 
            Chain = Tuple{ntuple(i -> Acyl, length(T))...}
        else
            n = startswith(rad, "d") ? 2 : 
            startswith(rad, "t") ? 3 :
            startswith(rad, "e") ? 4 : 1
            Chain = [Acyl for i in 1:length(T)]
            rrad = endswith(rad, "P-") ? Alkenyl : endswith(rad, "O-") ? Alkyl : Acyl
            for i in 1:n
                Chain[i] = rrad
            end
            Chain = Tuple{Chain...}
        end
        CarbonChain{Chain}(parse_singlechain(Radyl, cchain, pos, sil, mod)...)
    end
end

function parse_carbonchain(::Type{T}, cchain, pos, sil, mod) where {T <: Radyl}
    rad, cchain = match(r"([d,t,e]?[P,O]?-?)(\d+:\d+.*)", cchain)
    Chain = isnothing(rad) ? Acyl : isempty(rad) ? Acyl : endswith(rad, "P-") ? Alkenyl : endswith(rad, "O-") ? Alkyl : throw(ArgumentError("Invalid fattyacyl chain"))
    Chain <: T || throw(ArgumentError("Fattyacyl chain does not match to class"))
    CarbonChain{Chain}(parse_singlechain(Chain, cchain, pos, sil, mod)...)
end

function parse_carbonchain(::Type{T}, cchain, pos, sil, mod) where {T <: SPB}
    rad, cchain = match(r"([d,t,e]?[P,O]?-?)(\d+:\d+.*)", cchain)
    isempty(rad) || throw(ArgumentError("Invalid fattyacyl chain"))
    CarbonChain{T}(parse_singlechain(T, cchain, pos, sil, mod)...)
end

function parse_singlechain(::Type{T}, cchain, pos, sil, mod) where {T <: Radyl}
    cb, db = match(r"(\d+):(\d+)", cchain)
    cb = parse(UInt8, cb)
    db = parse(UInt8, db)
    isnothing(pos) && return cb, db, nothing, nothing
    ps = collect(eachmatch(r"(\d+)([EZ])?", pos))
    db = zeros(UInt8, length(ps))
    for (i, p) in enumerate(ps)
        x, e = p
        x = parse(UInt8, x) * 0x03
        e = isnothing(e) ? 0x00 : e == "Z" ? 0x01 : e == "E" ? 0x02 : throw(ArgumentError("Invalid fattyacyl chain"))
        db[i] = x + e
    end
    return cb, db, nothing, nothing
end

function parse_singlechain(::Type{T}, cchain, pos, sil, mod) where {T <: SPB}
    cb, db = match(r"(\d+):(\d+)", cchain)
    cb = parse(UInt8, cb)
    db = parse(UInt8, db)
    m = match(r";O(\d*)", mod)
    if isnothing(m)
        m = eachmatch(r"(\d+)OH", mod)
        if isempty(m)
            m = match(r";\(+OH\)+(\d*)", mod)
            n, = m
            sub = [Hydroxy() => isnothing(n) ? 0x01 : parse(UInt8, n)]
        else
            sub = [parse(UInt8, first(x.captures)) => Hydroxy() for x in m]
        end
    else
        n, = m
        sub = isnothing(n) ? 0x01 : parse(UInt8, n)
    end
    isnothing(pos) && return cb, db, sub, nothing
    ps = collect(eachmatch(r"(\d+)([EZ])?", pos))
    db = zeros(UInt8, length(ps))
    for (i, p) in enumerate(ps)
        x, e = p
        x = parse(UInt8, x) * 0x03
        e = isnothing(e) ? 0x00 : e == "Z" ? 0x01 : e == "E" ? 0x02 : throw(ArgumentError("Invalid fattyacyl chain"))
        db[i] = x + e
    end
    return cb, db, sub, nothing
end

function parse_sil(s)
    s = replace(s, r"(\d+[A-Z][a-z]*)" => s"[\1]")
    s = replace(s, r"(\d*),([\[,A-Z])" => s"\1\2")
end # TO FORMULA

function chemicalname(lipid::T) where {T <: Glycerophospholipid}
    position = interpretate_sn(lipid)
    sep = any(==(0), position) ? "_" : "/"
    string(class_abbr(lipid), " ", join(repr_singlechain.(lipid.chain), sep))
end

function chemicalname(lipid::T) where {T <: Sphingolipid}
    position = interpretate_position(lipid)
    sep = any(==(0), position) ? "_" : "/"
    string(class_abbr(lipid), " ", join(repr_singlechain.(lipid.chain), sep))
end

function interpretate_sn(lipid::T) where {T <: Lipid}
    position = zeros(Int, length(lipid.chain))
    sn = lipid.sn
    ep = length(lipid.chain) - 1
    bs = length(snposition(T)) + 1
    for i in eachindex(position)
        p, sn = divrem(sn, bs ^ ep)
        position[i] = p
        ep -= 1
    end
    position
end

function interpretate_position(lipid::T) where {T <: Sphingolipid}
    position = zeros(Int, length(lipid.chain))
    sn = lipid.position
    ep = length(lipid.chain) - 1
    bs = length(snposition(T)) + 1
    for i in eachindex(position)
        p, sn = divrem(sn, bs ^ ep)
        position[i] = p
        ep -= 1
    end
    position
end

function interpretate_db(db::UInt8)
    a, b = divrem(db, 3)
    string(a, b == 0 ? "" : b == 1 ? "Z" : b == 2 ? "E" : throw(ArgumentError("Invalid fattyacyl chain")))
end

function interpretate_sub(sub::UInt8)
    sub == 0 ? "" : string(";O", sub > 1 ? Int(sub) : "")
end 

function interpretate_sub(sub::Vector{<: Pair{UInt8, <: FunctionalGroup}})
    isempty(sub) && return ""
    sub = sort(sub; by = sub_abbr ∘ last)
    s = ""
    prev = ""
    for x in sub
        p, next = x
        next = sub_abbr(next)
        if next == prev
            s *= string(",", Int(p), next)
        else
            s *= string(";", Int(p), next)
            prev = next
        end
    end
    s
end 

class_abbr(::Ceramide) = "Cer"
class_abbr(::Phosphatidylcholine) = "PC"
class_abbr(::Phosphatidylethanolamine) = "PE"
class_abbr(::Cardiolipin) = "CL"
sub_abbr(::Hydroxy) = "OH"
repr_singlechain(c::CarbonChain{<: SPB}) = string(ncarbon(c), ":", ndoublebond(c), c.doublebond isa Vector ? string("(", join([interpretate_db(x) for x in c.doublebond], ","), ")") : "", interpretate_sub(c.substituent)) 
repr_singlechain(c::CarbonChain{<: Acyl}) = string(ncarbon(c), ":", ndoublebond(c), c.doublebond isa Vector ? string("(", join([interpretate_db(x) for x in c.doublebond], ","), ")") : "") 
repr_singlechain(c::CarbonChain{<: Alkyl}) = string("O-", ncarbon(c), ":", ndoublebond(c), c.doublebond isa Vector ? string("(", join([interpretate_db(x) for x in c.doublebond], ","), ")") : "") 
repr_singlechain(c::CarbonChain{<: Alkenyl}) = string("P-", ncarbon(c), ":", ndoublebond(c), c.doublebond isa Vector ? string("(", join([interpretate_db(x) for x in c.doublebond], ","), ")") : "") 

function repr_smiles(lipid::Glycerophospholipid)
    position = interpretate_sn(lipid)
    p = sortperm(position)
    cs = [string("(", repr_smiles(c), ")") for c in lipid.chain[p]]
    smi = repr_smiles(lipid.backbone)
    r = collect(eachmatch(r"\(O\)", smi))
    id = r[end - length(cs) + 1].match.offset
    s = smi[begin:id]
    prev = -3
    for (x, y) in zip(r[end - length(cs) + 1:end], cs)
        next = x.match.offset
        s *= smi[prev + 4:next]
        prev = next
        s *= y
    end
    s
end

function repr_smiles(back::DehydratedChemical)
    "OC(O)C(O)"
end
repr_smiles(chain::CarbonChain{Acyl, UInt8}) = "OC(=O)" * "C" ^ (ncarbon(chain) - 1)
function repr_smiles(chain::CarbonChain{Acyl})
    i = 2
    pos = sort!([divrem(db, 3) for db in chain.doublebond]; by = first)
    j = 1
    s = "OC(=O)"
    while i <= ncarbon(chain)
        if j <= lastindex(pos) && i == first(pos[j])
            if last(pos[j]) == 0
                s *= "C=C"
            elseif last(pos[j]) == 1
                if last(s) == '\\'
                    s *= "C=C/"
                elseif last(s) == '/'
                    s *= "C=C\\"
                else
                    s *= "\\C=C/"
                end
            elseif last(pos[j]) == 2
                if last(s) == '\\'
                    s *= "C=C\\"
                elseif last(s) == '/'
                    s *= "C=C/"
                else
                    s *= "\\C=C\\"
                end
            end
            j += 1
            i += 2
        else
            s *= "C"
            i += 1
        end
    end
    s
end

function repr_smiles(chain::CarbonChain{Alkyl, UInt8})
    "O" * "C" ^ ncarbon(chain)
end
function repr_smiles(chain::CarbonChain{Alkyl})
    i = 2
    pos = sort!([divrem(db, 3) for db in chain.doublebond]; by = first)
    j = 1
    s = "OC"
    while i <= ncarbon(chain)
        if j <= lastindex(pos) && i == first(pos[j])
            if last(pos[j]) == 0
                s *= "C=C"
            elseif last(pos[j]) == 1
                if last(s) == '\\'
                    s *= "C=C/"
                elseif last(s) == '/'
                    s *= "C=C\\"
                else
                    s *= "\\C=C/"
                end
            elseif last(pos[j]) == 2
                if last(s) == '\\'
                    s *= "C=C\\"
                elseif last(s) == '/'
                    s *= "C=C/"
                else
                    s *= "\\C=C\\"
                end
            end
            j += 1
            i += 2
        else
            s *= "C"
            i += 1
        end
    end
    s
end

function repr_smiles(chain::CarbonChain{Alkenyl, UInt8})
    "OC=C" * "C" ^ (ncarbon(chain) - 1)
end
function repr_smiles(chain::CarbonChain{Alkenyl})
    i = 3
    pos = sort!([divrem(db, 3) for db in chain.doublebond]; by = first)
    j = 1
    s = "O\\C=C/"
    while i <= ncarbon(chain)
        if j <= lastindex(pos) && i == first(pos[j])
            if last(pos[j]) == 0
                s *= "C=C"
            elseif last(pos[j]) == 1
                if last(s) == '\\'
                    s *= "C=C/"
                elseif last(s) == '/'
                    s *= "C=C\\"
                else
                    s *= "\\C=C/"
                end
            elseif last(pos[j]) == 2
                if last(s) == '\\'
                    s *= "C=C\\"
                elseif last(s) == '/'
                    s *= "C=C/"
                else
                    s *= "\\C=C\\"
                end
            end
            j += 1
            i += 2
        else
            s *= "C"
            i += 1
        end
    end
    s
end