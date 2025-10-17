using CellMetabolismBase

glycolysis_pathway = MetabolicPathway(
    (:Glucose_media, :Lactate_media),
    (
        (:GLUT1, (:Glucose_media,), (:Glucose,)),
        (:HK1, (:Glucose, :ATP), (:G6P, :ADP), (:Phosphate,), (:G6P,)),
        (:GPI, (:G6P,), (:F6P,)),
        (:PFKP, (:F6P, :ATP), (:F16BP, :ADP), (:Phosphate, :ADP, :F26BP), (:ATP, :Citrate)),
        (:ALDO, (:F16BP,), (:GAP, :DHAP)),
        (:TPI, (:DHAP,), (:GAP,)),
        (:GAPDH, (:GAP, :NAD, :Phosphate), (:BPG, :NADH)),
        (:PGK, (:BPG, :ADP), (:ThreePG, :ATP)),
        (:PGAM, (:ThreePG,), (:TwoPG,)),
        (:ENO, (:TwoPG,), (:PEP,)),
        (:PKM2, (:PEP, :ADP), (:Pyruvate, :ATP), (:F16BP,), (:Phenylalanine,)),
        (:LDH, (:Pyruvate, :NADH), (:Lactate, :NAD)),
        (:MCT, (:Lactate,), (:Lactate_media,)),
        (:AK, (:ADP, :ADP), (:ATP, :AMP)),
        (:NDPK, (:ATP, :NDP), (:ADP, :NTP)),
        (:CK, (:ATP, :Creatine), (:Phosphocreatine, :ADP)),
        (:ATPase, (:ATP,), (:ADP, :Phosphate)),
    ),
)
