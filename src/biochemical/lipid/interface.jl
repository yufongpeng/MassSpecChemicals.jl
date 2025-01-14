originalmolecule(chain::CarbonChain{Alkyl}) = FattyAlcohol(chain)
originalmolecule(chain::CarbonChain{Alkenyl}) = FattyAlcohol(chain)
originalmolecule(chain::CarbonChain{Acyl}) = FattyAcid(chain)
originalmolecule(chain::CarbonChain{SPB}) = SphingoBone(nothing, chain, 0x00)
leavinggroup(::CarbonChain{Acyl}) = Dehydroxy()
leavinggroup(::CarbonChain{SPB}) = Dehydrogen()
originalmolecule(chain::CarbonChain{<: AbstractSTRing}) = SterolBone(chain)
leavinggroup(::CarbonChain{<: AbstractSTRing}) = Dehydrogen()

# interface MassSpecChemicals.BioChemicals.Lipids for new lipid subclass