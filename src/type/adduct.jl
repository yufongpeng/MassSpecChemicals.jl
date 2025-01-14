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
"""
struct PosAdduct <: AbstractPosAdduct
    kmer::Int
    formula::String
    charge::Int
end

"""
    NegAdduct <: AbstractNegAdduct

Adduct with negative charges.
"""
struct NegAdduct <: AbstractNegAdduct
    kmer::Int
    formula::String
    charge::Int
end

struct LossElectron <: AbstractPosAdduct end
struct Protonation <: AbstractPosAdduct end
struct ProtonationNLH2O <: AbstractPosAdduct end
struct ProtonationNL2H2O <: AbstractPosAdduct end
struct ProtonationNL3H2O <: AbstractPosAdduct end
struct DiProtonation <: AbstractPosAdduct end
struct TriProtonation <: AbstractPosAdduct end
struct AddNH4 <: AbstractPosAdduct end
struct AddHNH4 <: AbstractPosAdduct end
struct Add2NH4 <: AbstractPosAdduct end
struct Sodization <: AbstractPosAdduct end
struct SodizationProtonation <: AbstractPosAdduct end
struct DiSodization <: AbstractPosAdduct end
struct SodizationAddNH4 <: AbstractPosAdduct end
struct Potassiation <: AbstractPosAdduct end
struct PotassiationProtonation <: AbstractPosAdduct end
struct DiPotassiation <: AbstractPosAdduct end
struct Lithiation <: AbstractPosAdduct end
struct LithiationProtonation <: AbstractPosAdduct end
struct DiLithiation <: AbstractPosAdduct end
struct Silveration <: AbstractPosAdduct end

struct AddElectron <: AbstractNegAdduct end
struct Deprotonation <: AbstractNegAdduct end
struct DeprotonationNLH2O <: AbstractNegAdduct end
struct DiDeprotonation <: AbstractNegAdduct end
struct TriDeprotonation <: AbstractNegAdduct end
struct AddOAc <: AbstractNegAdduct end
struct AddOFo <: AbstractNegAdduct end

struct LossCH2O <: AbstractNegAdduct end
struct AddO <: AbstractNegAdduct end
struct AddC2H2O <: AbstractNegAdduct end
struct LossCH8NO <: AbstractNegAdduct end
struct LossC2H8NO <: AbstractNegAdduct end
struct AddC3H5NO <: AbstractNegAdduct end
struct AddC2H5NO <: AbstractNegAdduct end
struct LossCH3 <: AbstractNegAdduct end
struct DeprotonationLossSerineAddH2O <: AbstractNegAdduct end
struct AddFluorine <: AbstractNegAdduct end
struct AddChlorine <: AbstractNegAdduct end
