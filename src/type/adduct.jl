"""
    AbstractAdduct

Abstract type for all kinds of adducts.
"""
abstract type AbstractAdduct end
"""
    AbstractPosAdduct <: AbstractAdduct

Abstract type for adducts with positive charges.
"""
abstract type AbstractPosAdduct <: AbstractAdduct end
"""
    AbstractNegAdduct <: AbstractAdduct

Abstract type for adducts with negative charges.
"""
abstract type AbstractNegAdduct <: AbstractAdduct end

"""
    PosAdduct <: AbstractPosAdduct

Adduct with positive charges. 

This general adduct cannot be applied to isotopic labeled chemicals which the adduct ion formed with isotopes loss. Adducts with isotopic labeled have to be defined as different adduct.

It is recommended to define custumized adduct when this condition may occur for more generalized use.
"""
struct PosAdduct <: AbstractPosAdduct
    kmer::Int
    formula::String
    ncharge::Int
end

"""
    NegAdduct <: AbstractNegAdduct

Adduct with negative charges. 

This general adduct cannot be applied to isotopic labeled chemicals which the adduct ion formed with isotopes loss. Adducts with isotopic labeled have to be defined as different adduct.

It is recommended to define custumized adduct when this condition may occur for more generalized use.
"""
struct NegAdduct <: AbstractNegAdduct
    kmer::Int
    formula::String
    ncharge::Int
end

"""
    LossElectron <: AbstractPosAdduct

[M]+
"""
struct LossElectron <: AbstractPosAdduct end
"""
    Protonation <: AbstractPosAdduct

[M+H]+
"""
struct Protonation <: AbstractPosAdduct end
"""
    ProtonationNLH2O <: AbstractPosAdduct

[M+H-H2O]+
"""
struct ProtonationNLH2O <: AbstractPosAdduct end
"""
    ProtonationNL2H2O <: AbstractPosAdduct

[M+H-2H2O]+
"""
struct ProtonationNL2H2O <: AbstractPosAdduct end
"""
    ProtonationNL3H2O <: AbstractPosAdduct

[M+H-3H2O]+
"""
struct ProtonationNL3H2O <: AbstractPosAdduct end
"""
    DiProtonation <: AbstractPosAdduct

[M+2H+]2+
"""
struct DiProtonation <: AbstractPosAdduct end
"""
    TriProtonation <: AbstractPosAdduct

[M+3H]3+
"""
struct TriProtonation <: AbstractPosAdduct end
"""
    AddNH4 <: AbstractPosAdduct

[M+NH4]+
"""
struct AddNH4 <: AbstractPosAdduct end
"""
    AddHNH4 <: AbstractPosAdduct

[M+H+NH4]2+
"""
struct AddHNH4 <: AbstractPosAdduct end
"""
    Add2NH4 <: AbstractPosAdduct

[M+2NH4]2+
"""
struct Add2NH4 <: AbstractPosAdduct end
"""
    Sodization <: AbstractPosAdduct

[M+Na]+
"""
struct Sodization <: AbstractPosAdduct end
"""
    SodizationProtonation <: AbstractPosAdduct

[M+Na+H]2+
"""
struct SodizationProtonation <: AbstractPosAdduct end
"""
    DiSodization <: AbstractPosAdduct

[M+2Na]2+
"""
struct DiSodization <: AbstractPosAdduct end
"""
    SodizationAddNH4 <: AbstractPosAdduct

[M+Na+NH4]2+
"""
struct SodizationAddNH4 <: AbstractPosAdduct end
"""
    Potassiation <: AbstractPosAdduct

[M+K]+
"""
struct Potassiation <: AbstractPosAdduct end
"""
    PotassiationProtonation <: AbstractPosAdduct

[M+K+H]2+
"""
struct PotassiationProtonation <: AbstractPosAdduct end
"""
    DiPotassiation <: AbstractPosAdduct

[M+2K]2+
"""
struct DiPotassiation <: AbstractPosAdduct end
"""
    Lithiation <: AbstractPosAdduct

[M+Li]+
"""
struct Lithiation <: AbstractPosAdduct end
"""
    LithiationProtonation <: AbstractPosAdduct

[M+Li+H]2+
"""
struct LithiationProtonation <: AbstractPosAdduct end
"""
    DiLithiation <: AbstractPosAdduct

[M+2Li]2+
"""
struct DiLithiation <: AbstractPosAdduct end
"""
    Silveration <: AbstractPosAdduct

[M+Ag]+
"""
struct Silveration <: AbstractPosAdduct end

"""
    AddElectron <: AbstractNegAdduct

[M]-
"""
struct AddElectron <: AbstractNegAdduct end
"""
    Deprotonation <: AbstractNegAdduct

[M-H]-
"""
struct Deprotonation <: AbstractNegAdduct end
"""
    DeprotonationNLH2O <: AbstractNegAdduct

[M-H-H2O]-
"""
struct DeprotonationNLH2O <: AbstractNegAdduct end
"""
    DiDeprotonation <: AbstractNegAdduct

[M-2H]2-
"""
struct DiDeprotonation <: AbstractNegAdduct end
"""
    TriDeprotonation <: AbstractNegAdduct

[M-3H]3-
"""
struct TriDeprotonation <: AbstractNegAdduct end
"""
    AddOAc <: AbstractNegAdduct

[M+OAc]-
"""
struct AddOAc <: AbstractNegAdduct end
"""
    AddOFo <: AbstractNegAdduct

[M+OFo]-
"""
struct AddOFo <: AbstractNegAdduct end

# struct LossCH2O <: AbstractNegAdduct end
# struct AddO <: AbstractNegAdduct end
# struct AddC2H2O <: AbstractNegAdduct end
# struct LossCH8NO <: AbstractNegAdduct end
# struct LossC2H8NO <: AbstractNegAdduct end
# struct AddC3H5NO <: AbstractNegAdduct end
# struct AddC2H5NO <: AbstractNegAdduct end
# struct LossCH3 <: AbstractNegAdduct end
# struct DeprotonationLossSerineAddH2O <: AbstractNegAdduct end
"""
    Fluoridation <: AbstractNegAdduct

[M+F]-
"""
struct Fluoridation <: AbstractNegAdduct end
"""
    Chloridation <: AbstractNegAdduct

[M+Cl]-
"""
struct Chloridation <: AbstractNegAdduct end
"""
    Demethylation <: AbstractNegAdduct

[M-Me]-
"""
struct Demethylation <: AbstractNegAdduct end
