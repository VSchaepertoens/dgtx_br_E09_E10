% https://github.com/neel-lab/DrawGlycan-SNFGv2
% cd to DrawGlycan-SNFGv2-master

%% draw individual N-glycans
drawglycan('Gal(b4)GlcNAc(b2)Man(a6)[Gal(b4 -U "3S")GlcNAc(b2)Man(a3)]Man(b4)GlcNAc(b4)[Fuc(a6)]GlcNAc(b)Asn(-CHAR)')
drawglycan('{Gal(-U "1 x")}{GlcNAc(b2)Man(a6)[GlcNAc(b2)Man(a3)]Man(b4)GlcNAc(b4)[Fuc(a6)]GlcNAc(b)Asn(-CHAR)}')
drawglycan('GlcNAc(b2)Man(a6)[GlcNAc(b2)Man(a3)]Man(b4)GlcNAc(b4)[Fuc(a6)]GlcNAc(b)Asn(-CHAR)')
drawglycan('GlcNAc(b2)Man(a6)[GlcNAc(b2)Man(a3)]Man(b4)GlcNAc(b4)GlcNAc(b)Asn(-CHAR)')

drawglycan('Man(a6)[GlcNAc(b2)Man(a3)]Man(b4)GlcNAc(b4)[Fuc(a6)]GlcNAc(b)Asn(-CHAR)')
drawglycan('Man(a6)[GlcNAc(b2)Man(a3)]Man(b4)GlcNAc(b4)GlcNAc(b)Asn(-CHAR)')

%% draw multiple glycans at the same time
drawglycanTile({'Gal(b4)GlcNAc(b2)Man(a6)[Gal(b4 -U "3S")GlcNAc(b2)Man(a3)]Man(b4)GlcNAc(b4)[Fuc(a6)]GlcNAc(b)Asn(-CHAR)', ...
    '{Gal(-U "1 x")}{GlcNAc(b2)Man(a6)[GlcNAc(b2)Man(a3)]Man(b4)GlcNAc(b4)[Fuc(a6)]GlcNAc(b)Asn(-CHAR)}'})
drawglycanTile({'GlcNAc(b2)Man(a6)[GlcNAc(b2)Man(a3)]Man(b4)GlcNAc(b4)[Fuc(a6)]GlcNAc(b)Asn(-CHAR)', ...
    'Man(a6)[GlcNAc(b2)Man(a3)]Man(b4)GlcNAc(b4)[Fuc(a6)]GlcNAc(b)Asn(-CHAR)'})