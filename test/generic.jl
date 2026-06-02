cglc = Chemical("Glucose", ["C" => 6, "H" => 12, "O" => 6]; retentiontime = 1.5, abbreviation = "Glc", SMILES = "")
cgld = Chemical("Glucose_d6", "C6H6D6O6"; retentiontime = 1.5, abbreviation = "Glc[D6]", SMILES = "")
cps = Chemical("PS 18:0/20:4(5Z,8Z,11Z,14Z)", "C44H78NO10P"; retentiontime = 7.8)
cpsi1 = Chemical("PS[D3,13C3] 18:0/20:4(5Z,8Z,11Z,14Z)", "C41[13C]3H75D3NO10P"; retentiontime = 7.8)
cpsi2 = Chemical("PS 18:0[D5]/20:4(5Z,8Z,11Z,14Z)", "C44H73D5NO10P"; retentiontime = 7.8)
fa1 =  Chemical("FA 18:0", "C18H36O2")
fa1i =  Chemical("FA 18:0[D5]", "C13D5H36O2")
cserine = Chemical("Serine", "C3H5NO2"; abbreviation = "Ser")
cserinei = Chemical("Serine[D3,13C3]", "[13C]3H2D3NO2"; abbreviation = "Ser[D3,13C3]")
lossserine = ChemicalLoss(cserine)
lossserinei = ChemicalLoss(cserinei)
losshserine = ChemicalLoss(AdductIon(cserine, "[M+H]+"))
losshserinei = ChemicalLoss(AdductIon(cserinei, "[M+H]+"))

# Generic structure interface
push!(cps.property, :schema => [
    ChemicalLoss(Proton()) => [
        lossserine => losshserine
    ] 
])
push!(cpsi1.property, :structure => [
    nothing => [
        lossserine => lossserinei,
        losshserine => losshserinei
    ],
    ChemicalLoss(Proton()) => [
        lossserine => lossserinei  
    ]
    ]
)
push!(cpsi1.property, :schema => [
    ChemicalLoss(Proton()) => [
        lossserine => losshserine
    ] 
])
push!(fa1.property, :chemicalscheme => [
    ChemicalLoss(Proton()) => :sn1fa
])
push!(fa1i.property, :chemicalscheme => [
    ChemicalLoss(Proton()) => :sn1fa
])
push!(cpsi2.property, :structure => [
    ChemicalLoss(Proton()) => [
        :sn1fa => AdductIon(fa1i, ChemicalLoss(Proton()))     
    ],
    losshserine => [
        :sn1fa => AdductIon(fa1i, ChemicalLoss(Proton()))     
    ]
])

# dimh = Adduct(2, "+H", 1)
# test all default adduct 
# icglcall = [AdductIon(cglc, k) for k in keys(MSC.ADDUCT_NAME)]

icglc = [AdductIon(cglc, "[M+H]+"), AdductIon(cglc, "[M+H-H2O]+")]
icgld = [parse_chemical(ChemicalEntityParser(ChemicalParser()), "[Glucose_d6+H]+"; formula = "C6H6D6O6", retentiontime = 1.5), AdductIon(parse_chemical(ChemicalParser(; formula = "C6H6D6O6"), "Glucose_d6"; retentiontime = 1.5), "[M+H]+"), AdductIon(cgld, "[M+H]+")]
icps = [ionize(cps; adduct = losshserine), AdductIon(cps, "[2M+H]+"), AdductIon(cps, "[M-H]-")]
icpsi1 = [ionize(cpsi1; adduct = losshserine), AdductIon(cpsi1, "[2M+H]+"), AdductIon(cpsi1, "[M-H]-")]
icpsi2 = [ionize(cpsi2; adduct = losshserine), AdductIon(cpsi2, "[2M+H]+"), AdductIon(cpsi2, "[M-H]-")]
cp1 = ChemicalSeries(icps[1], AdductIon(fa1, ChemicalLoss(Proton())))
cp2 = ChemicalSeries(icps[3] => lossserine => outputchemical(cp1))
cp3 = ChemicalSeries(icps[3], cp1)
cp4 = ChemicalSeries(isotopomerize(icps[3], ["[13C]" => 5]) => isotopomerize(completescheme(icps[1], outputchemical(cp1)), ["[13C]" => 5]))
cp5 = ChemicalSeries(isotopomerize(icpsi1[3], ["[13C]" => 5, "D" => 2]), ChemicalSeries(isotopomerize(completescheme(icpsi1[3], lossserine), ["D" => 1]), isotopomerize(completescheme(icpsi1[1], outputchemical(cp1)), ["[13C]" => 5])))
pt1 = peak_table(MSScan(Isotopologues(icglc[1]; abundance = 1e5, threshold = crit(1e1, 1e-2))))
