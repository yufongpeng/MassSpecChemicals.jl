exp1 = Table(;
    Chemical = AdductIon.([
        Chemical("Fructose", "C6H12O6"; retentiontime = 1.25),
        Chemical("Glucose", "C6H12O6"; retentiontime = 1.51), 
        Chemical("Galactose", "C6H12O6"; retentiontime = 2.005), 
        Chemical("Mannose", "C6H12O6"; retentiontime = 1.97), 
        Chemical("Glucosamine", "C6H13O5N"; retentiontime = 2.12),
        Chemical("Galactosamine", "C6H13O5N"; retentiontime = 2.61),
        Chemical("Mannosamine", "C6H13O5N"; retentiontime = 2.51),
    ], "[M+H]+"),
    Abundance1 = [100, 10000, 5000, 4000, 2000, 1000, 700], 
)

ci1 = CoelutingIsobars(
        [retentiontime => acrit(0.15)], 
        [
            TOF([]; resolution = 20000) => rcrit(1e-5)
        ], 
        exp1
    )
ib1 = isobar_table(ci1; threshold = rcrit(1e-5))

FA180 = AdductIon(Chemical("FA 18:0", "C18H36O2"), "[M-H]-")
FA181 = AdductIon(Chemical("FA 18:1", "C18H34O2"), "[M-H]-")
FA200 = AdductIon(Chemical("FA 20:0", "C20H40O2"), "[M-H]-")
exp2 = Table(;
    Chemical = ChemicalSeries.([
        AdductIon(Chemical("PC O-18:1/18:0", "C44H88NO7P"; retentiontime = 7.8), "[M+OFo]-") => FA180,
        AdductIon(Chemical("PC O-18:0/18:1", "C44H88NO7P"; retentiontime = 7.85), "[M+OFo]-") => FA181,
        AdductIon(Chemical("PS 18:0/20:0", "C44H86NO10P"; retentiontime = 7.7), "[M-H]-") => FA180,
        AdductIon(Chemical("PS 18:0/20:0", "C44H86NO10P"; retentiontime = 7.7), "[M-H]-") => FA200,
        AdductIon(Chemical("PC P-18:0/18:0", "C44H88NO7P"; retentiontime = 7.94), "[M+OFo]-") => FA180,
        AdductIon(Chemical("PC P-18:0/18:1", "C44H86NO7P"; retentiontime = 7.68), "[M+OFo]-") => FA181,
    ]),
    Abundance1 = [2000, 1500, 10000, 5000, 4000, 3000], 
)

ci2 = CoelutingIsobars(
        [retentiontime => acrit(0.15)], 
        [
            Quadrupole([]; accuracy = 0.1, fwhm = 1.3, digits = 1) => rcrit(1e-3),
            Quadrupole([]; accuracy = 0.1, fwhm = 1.3, digits = 1) => rcrit(1e-3)
        ], 
        exp2
    )
ib2 = isobar_table(ci2)