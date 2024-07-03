isadductequal(x::AbstractAdduct, y::AbstractAdduct) = isequal(x, y)

kmer(adduct::T) where {T <: AbstractAdduct} = hasfield(T, :kmer) ? adduct.kmer : 1
adductformula(adduct::T) where {T <: AbstractAdduct} = hasfield(T, :formula) ? adduct.formula : ""
charge(adduct::T) where {T <: AbstractAdduct} = hasfield(T, :charge) ? adduct.charge : 1
function adductelement(adduct::AbstractAdduct)
    el = Pair{SubString{String}, Int64}[]
    pos_adds = split(adductformula(adduct), "+", keepempty = false)
    for pos_add in pos_adds
        # M-H / M+FA-H
        neg_adds = split(pos_add, "-")
        pos_add = popfirst!(neg_adds)
        npos_add = match(r"^\d*", pos_add)
        npos_add = isempty(npos_add) ? 1 : parse(Int, npos_add.match)
        if npos_add > 1
            for e in parse_compound(string(pos_add))
                push!(el, first(e) => (last(e) * npos_add))
            end
        else
            for e in parse_compound(string(pos_add))
                push!(el, e)
            end
        end
        for neg_add in neg_adds
            nneg_add = match(r"^\d*", neg_add)
            nneg_add = isempty(nneg_add) ? 1 : parse(Int, nneg_add.match)
            for e in parse_compound(string(neg_add))
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
adductformula(::AddElectron) = ""
adductformula(::Deprotonation) = "-H"
adductformula(::DeprotonationNLH2O) = "-H-H2O"
adductformula(::DiDeprotonation) = "-2H"
adductformula(::TriDeprotonation) = "-3H"
adductformula(::AddOAc) = "+CH3COO"
adductformula(::AddHCOO) = "+HCOO"
adductformula(::LossCH2O) = "-CH2O"
adductformula(::AddO) = "+O"
adductformula(::AddC2H2O) = "+C2H2O"
adductformula(::LossCH8NO) = "-CH8NO"
adductformula(::LossC2H8NO) = "-C2H8NO"
adductformula(::AddC3H5NO) = "+C3H5NO"
adductformula(::AddC2H5NO) = "+C2H5NO"
adductformula(::LossCH3) = "-CH3"
adductformula(::DeprotonationLossSerineAddH2O) = "-C3H6NO2"

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
charge(::AddElectron) = 1
charge(::Deprotonation) = 1
charge(::DeprotonationNLH2O) = 1
charge(::DiDeprotonation) = 2
charge(::TriDeprotonation) = 3
charge(::AddOAc) = 1
charge(::AddHCOO) = 1
charge(::LossCH2O) = 1
charge(::AddO) = 1
charge(::AddC2H2O) = 1
charge(::LossCH8NO) = 1
charge(::LossC2H8NO) = 1
charge(::AddC3H5NO) = 1
charge(::AddC2H5NO) = 1
charge(::LossCH3) = 1
charge(::DeprotonationLossSerineAddH2O) = 1

adductelement(::LossElectron) = Pair{SubString{String}, Int64}[]
adductelement(::Protonation) = Pair{SubString{String}, Int64}["H" => 1]
adductelement(::ProtonationNLH2O) = Pair{SubString{String}, Int64}["H" => -1, "O" => -1]
adductelement(::ProtonationNL2H2O) = Pair{SubString{String}, Int64}["H" => -3, "O" => -2]
adductelement(::ProtonationNL3H2O) = Pair{SubString{String}, Int64}["H" => -5, "O" => -3]
adductelement(::DiProtonation) = Pair{SubString{String}, Int64}["H" => 2]
adductelement(::TriProtonation) = Pair{SubString{String}, Int64}["H" => 3]
adductelement(::AddNH4) = Pair{SubString{String}, Int64}["N" => 1, "H" => 4]
adductelement(::AddHNH4) = Pair{SubString{String}, Int64}["N" => 1, "H" => 5]
adductelement(::Add2NH4) = Pair{SubString{String}, Int64}["N" => 2, "H" => 8]
adductelement(::Sodization) = Pair{SubString{String}, Int64}["Na" => 1]
adductelement(::SodizationProtonation) = Pair{SubString{String}, Int64}["Na" => 1, "H" => 1]
adductelement(::DiSodization) = Pair{SubString{String}, Int64}["Na" => 2]
adductelement(::AddElectron) = Pair{SubString{String}, Int64}[]
adductelement(::Deprotonation) = Pair{SubString{String}, Int64}["H" => -1]
adductelement(::DeprotonationNLH2O) = Pair{SubString{String}, Int64}["H" => -3, "O" => -1]
adductelement(::DiDeprotonation) = Pair{SubString{String}, Int64}["H" => -2]
adductelement(::TriDeprotonation) = Pair{SubString{String}, Int64}["H" => -3]
adductelement(::AddOAc) = Pair{SubString{String}, Int64}["C" => 2, "H" => 3, "O" => 2]
adductelement(::AddHCOO) = Pair{SubString{String}, Int64}["C" => 1, "H" => 1, "O" => 2]
adductelement(::LossCH2O) = Pair{SubString{String}, Int64}["C" => -1, "H" => -2, "O" => -1]
adductelement(::AddO) = Pair{SubString{String}, Int64}["O" => 1]
adductelement(::AddC2H2O) = Pair{SubString{String}, Int64}["C" => 2, "H" => 2, "O" => 1]
adductelement(::LossCH8NO) = Pair{SubString{String}, Int64}["C" => -1, "H" => -8, "N" => -1, "O" => -1]
adductelement(::LossC2H8NO) = Pair{SubString{String}, Int64}["C" => -2, "H" => -8, "N" => -1, "O" => -1]
adductelement(::AddC3H5NO) = Pair{SubString{String}, Int64}["C" => 3, "H" => 5, "N" => 1, "O" => 1]
adductelement(::AddC2H5NO) = Pair{SubString{String}, Int64}["C" => 2, "H" => 5, "N" => 1, "O" => 1]
adductelement(::LossCH3) = Pair{SubString{String}, Int64}["C" => -1, "H" => -3]
adductelement(::DeprotonationLossSerineAddH2O) = Pair{SubString{String}, Int64}["C" => -3, "H" => -6, "N" => -1, "O" => -2]
