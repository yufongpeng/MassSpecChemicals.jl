cglc = Chemical("Glucose", ["C" => 6, "H" => 12, "O" => 6]; retentiontime = 1.5, abbreviation = "Glc", SMILES = "")
cgld = parse_chemical("Glucose-d6"; formula = "C6H6D6O6", retentiontime = 1.5, abbreviation = "Glc[D6]", SMILES = "")
cps = Chemical("PS 18:0/20:4(5Z,8Z,11Z,14Z)", "C44H78NO10P"; retentiontime = 7.8)
cpsi1 = Chemical("PS[D3,13C3] 18:0/20:4(5Z,8Z,11Z,14Z)", "C41[13C]3H75D3NO10P"; retentiontime = 7.8)
cpsi2 = Chemical("PS 18:0[D5]/20:4(5Z,8Z,11Z,14Z)", "C44H73D5NO10P"; retentiontime = 7.8)
lossserine = NegAdduct(1, "-C3H6NO2", 1)
cserine = Chemical("Serine", "C3H5NO2")
cserinei = Chemical("Serine[D3,13C3]", "[13C]3H2D3NO2")
clossserine = ChemicalLoss(cserine)

# lossserinei = NegAdduct(1, "-[13C]3H3D3NO2", 1)
push!(cpsi1.property, :adductisotopes => [lossserine => ["C" => 3, "[13C]" => -3, "H" => 3, "D" => -3]])
dimh = PosAdduct(2, "+H", 1)
# test all default adduct 
icglcall = [AdductIon(cglc, k) for k in keys(MSC.ADDUCT_NAME)]

icglc = [AdductIon(cglc, Protonation())]
icgld = [AdductIon("Glucose-d6", "[M+H]+"; formula = "C6H6D6O6", retentiontime = 1.5), AdductIon("Glucose-d6", Protonation(); formula = "C6H6D6O6", retentiontime = 1.5), AdductIon(cgld, Protonation())]
icps = [AdductIon(cps, lossserine), AdductIon(cps, dimh)]
icpsi1 = [AdductIon(cpsi1, lossserine), AdductIon(cpsi1, dimh)]
icpsi2 = [AdductIon(cpsi2, lossserine), AdductIon(cpsi2, dimh)]
cp1 = ChemicalPair(icps[1], AdductIon(Chemical("FA 20:4", "C20H32O2"; retentiontime = 7.78), Deprotonation()))
cp2 = ChemicalSeries(AdductIon(cps, "[M-H]-") => clossserine => cp1.product)
cp3 = ChemicalPair(AdductIon(cps, "[M-H]-"), cp1)
cp4 = ChemicalSeries(Isotopomers(AdductIon(cps, "[M-H]-"), ["[13C]" => 5]) => Isotopomers(cp1.product, ["[13C]" => 5]))
cp5 = ChemicalPair(Isotopomers(AdductIon(cpsi1, "[M-H]-"), ["[13C]" => 5, "D" => 2]), ChemicalPair(ChemicalLoss(Isotopomers(cserinei, ["D" => 1])), Isotopomers(cp1.product, ["[13C]" => 5])))
pt1 = peak_table(MSScan(Isotopologues(icglc[1]; abundance = 1e5, threshold = crit(1e1, 1e-2))))
