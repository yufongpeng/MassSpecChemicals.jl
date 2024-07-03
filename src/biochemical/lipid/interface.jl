originalmolecule(chain::CarbonChain{Alkyl}) = FattyAlcohol(chain)
originalmolecule(chain::CarbonChain{Alkenyl}) = FattyAlcohol(chain)
originalmolecule(chain::CarbonChain{Acyl}) = FattyAcid(chain)
originalmolecule(chain::CarbonChain{SPB}) = SphingoBone(nothing, chain, 0x00)
leavinggroup(::CarbonChain{Acyl}) = Dehydroxy()
leavinggroup(::CarbonChain{SPB}) = Dehydrogen()
originalmolecule(chain::CarbonChain{<: AbstractSTRing}) = SterolBone(chain)
leavinggroup(::CarbonChain{<: AbstractSTRing}) = Dehydrogen()

snposition(::Type{<: Glycerophospholipid}) = ["sn-1", "sn-2"]
snposition(::Type{<: Glycerolipid}) = ["sn-1", "sn-2", "sn-3"]
snposition(::Type{<: Bisradylglycerophosphoglycerol}) = ["sn-1", "sn-2", "sn-1'", "sn-2'"]
snposition(::Type{<: Bisphosphatidicacid}) = ["sn-2", "sn-3", "sn-2'", "sn-3'"]
ncarbon(chain::CarbonChain) = chain.carbon
ndoublebond(chain::CarbonChain{S, UInt8}) where S = Int(chain.doublebond)
ndoublebond(chain::CarbonChain{S}) where S = length(chain.doublebond)