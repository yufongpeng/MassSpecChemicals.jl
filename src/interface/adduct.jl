"""
    isadductequal(x::AbstractAdduct, y::AbstractAdduct)

Determine whether two adducts are chemically equivalent. By default, it transforms both chemicals by `isadductequaltransform` and compares them by `istransformedadductequal`.
"""
isadductequal(x::AbstractAdduct, y::AbstractAdduct) = istransformedadductequal(isadductequaltransform(x), isadductequaltransform(y))

"""
    isadductequaltransform(x::AbstractAdduct) 

Return an object for comparison with other adducts by `istransformedadductequal`. 
"""
isadductequaltransform(x::AbstractAdduct) = x

"""
    istransformedadductequal(x::AbstractAdduct, y::AbstractAdduct)
    istransformedadductequal(x::PosAdduct, y::PosAdduct)
    istransformedadductequal(x::NegAdduct, y::NegAdduct)

Determine whether two chemicals are chemically equivalent after applying `isadductequaltransform`. It defaults to `isequal`.
"""
istransformedadductequal(x::AbstractAdduct, y::AbstractAdduct) = isequal(x, y)
istransformedadductequal(x::PosAdduct, y::PosAdduct) = kmer(x) == kmer(y) && ncharge(x) == ncharge(y) && Set(mapreduce(a -> split(a, "-"), vcat, split(adductformula(x), "+"))) == Set(mapreduce(a -> split(a, "-"), vcat, split(adductformula(y), "+")))
istransformedadductequal(x::NegAdduct, y::NegAdduct) = kmer(x) == kmer(y) && ncharge(x) == ncharge(y) && Set(mapreduce(a -> split(a, "-"), vcat, split(adductformula(x), "+"))) == Set(mapreduce(a -> split(a, "-"), vcat, split(adductformula(y), "+")))

"""
    kmer(adduct)

The number of core chemical. For instance, 2 for [2M+H]+.
"""
kmer(adduct::T) where {T <: AbstractAdduct} = hasfield(T, :kmer) ? adduct.kmer : 1

"""
    adductformula(adduct)

The formula for adduct. For instance,  `"-H"` for [M-H]-, `"+OAc"` for [M+OAc]-.
"""
adductformula(adduct::T) where {T <: AbstractAdduct} = hasfield(T, :formula) ? adduct.formula : nothing

"""
    charge(adduct)

The charge of adduct (positive or negative). For instance, -1 for [M-H]-, 2 for [M+2H]2+. The default value for positive and negative adduct are 1 and -1.
"""
charge(adduct::T) where {T <: AbstractPosAdduct} = 1
charge(adduct::T) where {T <: AbstractNegAdduct} = -1
charge(adduct::PosAdduct) = adduct.ncharge
charge(adduct::NegAdduct) = -1 * adduct.ncharge

"""
    ncharge(adduct)

The number of charges of adduct. For instance, 1 for [M-H]-, 2 for [M+2H]2+.
"""
ncharge(adduct::T) where {T <: AbstractAdduct} = abs(charge(adduct))
ncharge(adduct::PosAdduct) = adduct.ncharge
ncharge(adduct::NegAdduct) = adduct.ncharge

"""
    adductelements(adduct)

The elements changed with adduct formation. For generic adduct, it uses `adductformula` to calculate elements. If non-element strings are used in `adductformula`, defining custumized `adductelements` is required. 
"""
function adductelements(adduct::AbstractAdduct)
    el = Pair{String, Int}[]
    pos_adds = split(adductformula(adduct), "+", keepempty = false)
    for pos_add in pos_adds
        # M-H / M+FA-H
        neg_adds = split(pos_add, "-")
        pos_add = popfirst!(neg_adds)
        npos_add = match(r"^\d*", pos_add)
        npos_add = isempty(npos_add) ? 1 : parse(Int, npos_add.match)
        if npos_add > 1
            for e in chemicalelements(string(pos_add))
                push!(el, first(e) => (last(e) * npos_add))
            end
        else
            for e in chemicalelements(string(pos_add))
                push!(el, e)
            end
        end
        for neg_add in neg_adds
            nneg_add = match(r"^\d*", neg_add)
            nneg_add = isempty(nneg_add) ? 1 : parse(Int, nneg_add.match)
            for e in chemicalelements(string(neg_add))
                push!(el, first(e) => -last(e) * nneg_add)
            end
        end
    end
    el
end

adductformula(::LossElectron) = ""
adductformula(::Protonation) = "+H"
adductformula(::ProtonationNLH2O) = "+H-H2O"
adductformula(::ProtonationNL2H2O) = "+H-2H2O"
adductformula(::ProtonationNL3H2O) = "+H-3H2O"
adductformula(::DiProtonation) = "+2H"
adductformula(::TriProtonation) = "+3H"
adductformula(::AddNH4) = "+NH4"
adductformula(::AddHNH4) = "+H+NH4"
adductformula(::Add2NH4) = "+2NH4"
adductformula(::Sodization) = "+Na"
adductformula(::SodizationProtonation) = "+Na+H"
adductformula(::DiSodization) = "+2Na"
adductformula(::SodizationAddNH4) = "+Na+NH4"
adductformula(::Potassiation) = "+K"
adductformula(::PotassiationProtonation) = "+K+H"
adductformula(::DiPotassiation) = "+2K"
adductformula(::Lithiation) = "+Li"
adductformula(::LithiationProtonation) = "+Li+H"
adductformula(::DiLithiation) = "+2Li"
adductformula(::Silveration) = "+Ag"
adductformula(::AddElectron) = ""
adductformula(::Deprotonation) = "-H"
adductformula(::DeprotonationNLH2O) = "-H-H2O"
adductformula(::DiDeprotonation) = "-2H"
adductformula(::TriDeprotonation) = "-3H"
adductformula(::AddOAc) = "+OAc"
adductformula(::AddOFo) = "+OFo"
# adductformula(::LossCH2O) = "-CH2O"
# adductformula(::AddO) = "+O"
# adductformula(::AddC2H2O) = "+C2H2O"
# adductformula(::LossCH8NO) = "-CH8NO"
# adductformula(::LossC2H8NO) = "-C2H8NO"
# adductformula(::AddC3H5NO) = "+C3H5NO"
# adductformula(::AddC2H5NO) = "+C2H5NO"
# adductformula(::LossCH3) = "-CH3"
# adductformula(::DeprotonationLossSerineAddH2O) = "-H-Serine+H2O"
adductformula(::Fluoridation) = "+F"
adductformula(::Chloridation) = "+Cl"
adductformula(::Demethylation) = "-Me"

charge(::LossElectron) = 1
charge(::Protonation) = 1
charge(::ProtonationNLH2O) = 1
charge(::ProtonationNL2H2O) = 1
charge(::ProtonationNL3H2O) = 1
charge(::DiProtonation) = 2
charge(::TriProtonation) = 3
charge(::AddNH4) = 1
charge(::AddHNH4) = 2
charge(::Add2NH4) = 2
charge(::Sodization) = 1
charge(::SodizationProtonation) = 2
charge(::DiSodization) = 2
charge(::SodizationAddNH4) = 2
charge(::Potassiation) = 1
charge(::PotassiationProtonation) = 2
charge(::DiPotassiation) = 2
charge(::Lithiation) = 1
charge(::LithiationProtonation) = 2
charge(::DiLithiation) = 2
charge(::Silveration) = 1
charge(::AddElectron) = -1
charge(::Deprotonation) = -1
charge(::DeprotonationNLH2O) = -1
charge(::DiDeprotonation) = -2
charge(::TriDeprotonation) = -3
charge(::AddOAc) = -1
charge(::AddOFo) = -1
# charge(::LossCH2O) = -1
# charge(::AddO) = -1
# charge(::AddC2H2O) = -1
# charge(::LossCH8NO) = -1
# charge(::LossC2H8NO) = -1
# charge(::AddC3H5NO) = -1
# charge(::AddC2H5NO) = -1
# charge(::LossCH3) = -1
# charge(::DeprotonationLossSerineAddH2O) = -1
charge(::Fluoridation) = -1
charge(::Chloridation) = -1
charge(::Demethylation) = -1

adductelements(::LossElectron) = Pair{String, Int64}[]
adductelements(::Protonation) = ["H" => 1]
adductelements(::ProtonationNLH2O) = ["H" => -1, "O" => -1]
adductelements(::ProtonationNL2H2O) = ["H" => -3, "O" => -2]
adductelements(::ProtonationNL3H2O) = ["H" => -5, "O" => -3]
adductelements(::DiProtonation) = ["H" => 2]
adductelements(::TriProtonation) = ["H" => 3]
adductelements(::AddNH4) = ["N" => 1, "H" => 4]
adductelements(::AddHNH4) = ["N" => 1, "H" => 5]
adductelements(::Add2NH4) = ["N" => 2, "H" => 8]
adductelements(::Sodization) = ["Na" => 1]
adductelements(::SodizationProtonation) = ["Na" => 1, "H" => 1]
adductelements(::DiSodization) = ["Na" => 2]
adductelements(::SodizationAddNH4) = ["Na" => 1, "N" => 1, "H" => 4]
adductelements(::Potassiation) = ["K" => 1]
adductelements(::PotassiationProtonation) = ["K" => 1, "H" => 1]
adductelements(::DiPotassiation) = ["K" => 2]
adductelements(::Lithiation) = ["Li" => 1]
adductelements(::LithiationProtonation) = ["Li" => 1, "H" => 1]
adductelements(::DiLithiation) = ["Li" => 2]
adductelements(::Silveration) = ["Ag" => 1]
adductelements(::AddElectron) = Pair{String, Int64}[]
adductelements(::Deprotonation) = ["H" => -1]
adductelements(::DeprotonationNLH2O) = ["H" => -3, "O" => -1]
adductelements(::DiDeprotonation) = ["H" => -2]
adductelements(::TriDeprotonation) = ["H" => -3]
adductelements(::AddOAc) = ["C" => 2, "H" => 3, "O" => 2]
adductelements(::AddOFo) = ["C" => 1, "H" => 1, "O" => 2]
# adductelements(::LossCH2O) = ["C" => -1, "H" => -2, "O" => -1]
# adductelements(::AddO) = ["O" => 1]
# adductelements(::AddC2H2O) = ["C" => 2, "H" => 2, "O" => 1]
# adductelements(::LossCH8NO) = ["C" => -1, "H" => -8, "N" => -1, "O" => -1]
# adductelements(::LossC2H8NO) = ["C" => -2, "H" => -8, "N" => -1, "O" => -1]
# adductelements(::AddC3H5NO) = ["C" => 3, "H" => 5, "N" => 1, "O" => 1]
# adductelements(::AddC2H5NO) = ["C" => 2, "H" => 5, "N" => 1, "O" => 1]
# adductelements(::LossCH3) = ["C" => -1, "H" => -3]
# adductelements(::DeprotonationLossSerineAddH2O) = ["C" => -3, "H" => -6, "N" => -1, "O" => -2]
adductelements(::Fluoridation) = ["F" => 1]
adductelements(::Chloridation) = ["Cl" => 1]
adductelements(::Demethylation) = ["C" => -1, "H" => -3]

Broadcast.broadcastable(x::AbstractAdduct) = Ref(x)