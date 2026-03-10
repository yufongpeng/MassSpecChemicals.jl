"""
    AbstractAdduct

Abstract type for all kinds of adducts.
"""
abstract type AbstractAdduct end

"""
    ComposedAdduct{T} <: AbstractAdduct

Adducts composed of multiple adducts and chemical gain or loss.

# Fields
* `adducts::Vector{T}`

Attribute `kmer` is determined by the first adduct, and `charge` is determined by the sum of charges of all adducts.
"""
struct ComposedAdduct{T} <: AbstractAdduct
    adducts::Vector{T}
end 

"""
    Adduct <: AbstractAdduct

Gemeric adduct. 

# Fields
* `kmer::Int`
* `formula::String`
* `charge::Int`

It is recommended to define custumized adduct when this condition may occur for more custumized use.
"""
struct Adduct <: AbstractAdduct
    kmer::Int
    formula::String
    charge::Int
end

"""
    LossElectron <: AbstractAdduct

[M]+
"""
struct LossElectron <: AbstractAdduct end
"""
    Protonation <: AbstractAdduct

[M+H]+
"""
struct Protonation <: AbstractAdduct end
"""
    ProtonationNLH2O <: AbstractAdduct

[M+H-H2O]+
"""
struct ProtonationNLH2O <: AbstractAdduct end
"""
    ProtonationNL2H2O <: AbstractAdduct

[M+H-2H2O]+
"""
struct ProtonationNL2H2O <: AbstractAdduct end
"""
    ProtonationNL3H2O <: AbstractAdduct

[M+H-3H2O]+
"""
struct ProtonationNL3H2O <: AbstractAdduct end
"""
    DiProtonation <: AbstractAdduct

[M+2H+]2+
"""
struct DiProtonation <: AbstractAdduct end
"""
    TriProtonation <: AbstractAdduct

[M+3H]3+
"""
struct TriProtonation <: AbstractAdduct end
"""
    AddNH4 <: AbstractAdduct

[M+NH4]+
"""
struct AddNH4 <: AbstractAdduct end
"""
    AddNH4Protonation <: AbstractAdduct

[M+H+NH4]2+
"""
struct AddNH4Protonation <: AbstractAdduct end
"""
    Add2NH4 <: AbstractAdduct

[M+2NH4]2+
"""
struct Add2NH4 <: AbstractAdduct end
"""
    Sodization <: AbstractAdduct

[M+Na]+
"""
struct Sodization <: AbstractAdduct end
"""
    SodizationProtonation <: AbstractAdduct

[M+Na+H]2+
"""
struct SodizationProtonation <: AbstractAdduct end
"""
    DiSodization <: AbstractAdduct

[M+2Na]2+
"""
struct DiSodization <: AbstractAdduct end
"""
    SodizationAddNH4 <: AbstractAdduct

[M+Na+NH4]2+
"""
struct SodizationAddNH4 <: AbstractAdduct end
"""
    Potassiation <: AbstractAdduct

[M+K]+
"""
struct Potassiation <: AbstractAdduct end
"""
    PotassiationProtonation <: AbstractAdduct

[M+K+H]2+
"""
struct PotassiationProtonation <: AbstractAdduct end
"""
    DiPotassiation <: AbstractAdduct

[M+2K]2+
"""
struct DiPotassiation <: AbstractAdduct end
"""
    Lithiation <: AbstractAdduct

[M+Li]+
"""
struct Lithiation <: AbstractAdduct end
"""
    LithiationProtonation <: AbstractAdduct

[M+Li+H]2+
"""
struct LithiationProtonation <: AbstractAdduct end
"""
    DiLithiation <: AbstractAdduct

[M+2Li]2+
"""
struct DiLithiation <: AbstractAdduct end
"""
    Silveration <: AbstractAdduct

[M+Ag]+
"""
struct Silveration <: AbstractAdduct end

"""
    AddElectron <: AbstractAdduct

[M]-
"""
struct AddElectron <: AbstractAdduct end
"""
    Deprotonation <: AbstractAdduct

[M-H]-
"""
struct Deprotonation <: AbstractAdduct end
"""
    DeprotonationNLH2O <: AbstractAdduct

[M-H-H2O]-
"""
struct DeprotonationNLH2O <: AbstractAdduct end
"""
    DiDeprotonation <: AbstractAdduct

[M-2H]2-
"""
struct DiDeprotonation <: AbstractAdduct end
"""
    TriDeprotonation <: AbstractAdduct

[M-3H]3-
"""
struct TriDeprotonation <: AbstractAdduct end
"""
    AddOAc <: AbstractAdduct

[M+OAc]-
"""
struct AddOAc <: AbstractAdduct end
"""
    AddOFo <: AbstractAdduct

[M+OFo]-
"""
struct AddOFo <: AbstractAdduct end
"""
    Fluoridation <: AbstractAdduct

[M+F]-
"""
struct Fluoridation <: AbstractAdduct end
"""
    Chloridation <: AbstractAdduct

[M+Cl]-
"""
struct Chloridation <: AbstractAdduct end
"""
    Demethylation <: AbstractAdduct

[M-Me]-
"""
struct Demethylation <: AbstractAdduct end
