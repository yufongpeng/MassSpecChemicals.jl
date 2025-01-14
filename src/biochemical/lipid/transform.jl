dissociate_headgroup(lipid::Glycerolipid) = lipid
function dissociate_headgroup(lipid::Union{<: Glycerophospholipid, <: Omodifiedradylglycerol})
    sn = decode_sn(lipid)
    n = ncarbonchain(lipid)
    cls = @match n begin
        1   => Monoradylglycerol
        2   => Diradylglycerol
        3   => Triradylglycerol
    end
    sn = sn .* ((nchainposition(cls) + 1) .^ ((n - 1):-1:0))
    cls(last(getchaincomponent(lipid.backbone)), lipid.chain, UInt8(sum(sn)))
end

function dissociate_headgroup(lipid::CeramideBone)
    position = decode_position(lipid)
    isempty(position) && return SphingoBone(nothing, lipid.chain, 0x00)
    spb, acyl = lipid.chain
    sub = deepcopy(spb.substituent)
    for p in position
        push!(sub, UInt8(p) => Hydroxy())
    end
    SphingoBone(nothing, (make_carbonchain(SPB, spb.carbon, spb.doublebond, sort!(sub; by = x -> ((sub_abbr ∘ last)(x), first(x))), spb.isotopiclabel), acyl), 0x00)
end
dissociate_headgroup(lipid::Ceramide) = lipid
function dissociate_headgroup(lipid::SphingoidBaseBone)
    position = decode_position(lipid)
    isempty(position) && return SphingoBone(nothing, lipid.chain, 0x00)
    sub = deepcopy(lipid.chain.substituent)
    for p in position
        push!(sub, UInt8(p) => Hydroxy())
    end
    SphingoBone(nothing, make_carbonchain(SPB, lipid.chain.carbon, lipid.chain.doublebond, sort!(sub; by = x -> ((sub_abbr ∘ last)(x), first(x))), lipid.chain.isotopiclabel), 0x00)
end
dissociate_headgroup(lipid::SphingoidBase) = lipid