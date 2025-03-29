using CellMetabolismBase

glycolysis_pathway = MetabolicPathway(
    (:Glucose_media, :Lactate_media),
    (
        (:GLUT1, (:Glucose_media,), (:Glucose,)),
        (:HK1, (:Glucose, :ATP), (:G6P, :ADP)),
        (:GPI, (:G6P,), (:F6P,)),
        (:PFKP, (:F6P, :ATP), (:F16BP, :ADP)),
        (:ALDO, (:F16BP,), (:GAP, :DHAP)),
        (:TPI, (:DHAP,), (:GAP,)),
        (:GAPDH, (:GAP, :NAD, :Phosphate), (:BPG, :NADH)),
        (:PGK, (:BPG, :ADP), (:ThreePG, :ATP)),
        (:PGAM, (:ThreePG,), (:TwoPG,)),
        (:ENO, (:TwoPG,), (:PEP,)),
        (:PKM2, (:PEP, :ADP), (:Pyruvate, :ATP)),
        (:LDH, (:Pyruvate, :NADH), (:Lactate, :NAD)),
        (:MCT, (:Lactate,), (:Lactate_media,)),
        (:AK, (:ADP, :ADP), (:ATP, :AMP)),
        (:NDPK, (:ATP, :NDP), (:ADP, :NTP)),
        (:CK, (:ATP, :Creatine), (:Phosphocreatine, :ADP)),
        (:ATPase, (:ATP,), (:ADP, :Phosphate)),
    ),
)
