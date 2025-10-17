using CellMetabolismBase, LabelledArrays, Measurements

"Fraction of intracellular volume occupied by water with macromolecules making up the rest"
water_fraction_cell_volume = 0.66

"Fraction of intracellular volume occupied by cytosol"
cytosol_fraction_cell_volume = 0.66

"""
Correction for cytosolic volume and its water content to better estimate concentration of
cytosolic molecules from total cellular concentration and vice versa
"""
cell_volume_correction = water_fraction_cell_volume * cytosol_fraction_cell_volume

"Intracellular protein density, mg protein per µl of cell volume"
cell_protein_density = 0.2


glycolysis_init_conc_w_uncertainty = LVector(
    Glucose_media = 25e-3,
    Glucose = (4.1e-3 ± 8e-4) / cell_volume_correction,
    G6P = (1.7e-4 ± 3e-5) / cell_volume_correction,
    F6P = (9.8e-5 ± 1.0e-5) / cell_volume_correction,
    F16BP = (5.5e-4 ± 1.0e-4) / cell_volume_correction,
    GAP = (9.6e-5 ± 2.0e-5) / cell_volume_correction,
    DHAP = (8.2e-4 ± 1.9e-4) / cell_volume_correction,
    BPG = (1.1e-6 ± 3e-7) / cell_volume_correction,
    ThreePG = (1.6e-3 ± 6e-4) / cell_volume_correction,
    TwoPG = (4.9e-4 ± 2e-4) / cell_volume_correction,
    PEP = (9.7e-5 ± 1.9e-5) / cell_volume_correction,
    Pyruvate = (6.1e-4 ± 2.3e-4) / cell_volume_correction,
    Lactate = (9.6e-3 ± 2.0e-3) / cell_volume_correction,
    Lactate_media = 2e-3,
    ATP = (3.6e-3 ± 5e-4) / cell_volume_correction,
    ADP = (7.3e-4 ± 1.5e-4) / cell_volume_correction,
    AMP = (1.6e-4 ± 6e-5) / cell_volume_correction,
    Phosphate = (4.0e-3 ± 9e-4) / cell_volume_correction,
    NTP = (0.00186 ± 0.00022) / cell_volume_correction,
    NDP = (0.000354 ± 0.000058) / cell_volume_correction,
    Phosphocreatine = 0.003 / cell_volume_correction,
    Creatine = 0.0003 / cell_volume_correction,
    NAD = (9.2e-4 ± 4.1e-4) / cell_volume_correction,
    NADH = (8.4e-5 ± 3.6e-5) / cell_volume_correction,
    F26BP = 0.0 / cell_volume_correction,
    Citrate = 0.0 / cell_volume_correction,
    Phenylalanine = 0.0 / cell_volume_correction,
)

glycolysis_init_conc = Measurements.value.(glycolysis_init_conc_w_uncertainty)


glycolysis_params_w_uncertainty = LVector(
    GLUT1_Km_Glucose = 20e-3 ± 4e-3,
    GLUT1_Conc = (cell_protein_density / cell_volume_correction) * (1.4e-4 ± 4e-5),
    GLUT1_Vmax = 1100.0 ± 300.0,
    GLUT1_Keq = 1.0 ± 0.0,
    GLUT1_MW = 54084.0 / 1000,
    HK1_K_Glucose = 4.9e-8 ± 0.4e-8,
    HK1_K_a_ATP = 9.3e-7 ± 1e-7,
    HK1_β_Glucose_ATP = 0.0016 ± 0.0002,
    HK1_K_G6P = 480e-6 ± 130e-6,
    HK1_K_a_ADP = 390e-6 ± 50e-6,
    HK1_K_a_G6P_cat = 33e-6 ± 8e-6,
    HK1_K_i_G6P_reg = 6.9e-6 ± 0.9e-6,
    HK1_K_a_Pi = 460e-6 ± 120e-6,
    HK1_Conc = (cell_protein_density / cell_volume_correction) * (6.0e-4 ± 1.2e-4),
    HK1_Vmax = 53.0 ± 5.0,
    HK1_Keq = 2700.0 ± 800.0,
    HK1_MW = 102486.0 / 1000,
    GPI_Km_G6P = 330e-6 ± 100e-6,
    GPI_Km_F6P = 70e-6 ± 20e-6,
    GPI_Conc = (cell_protein_density / cell_volume_correction) * (1.4e-3 ± 1e-4),
    GPI_Vmax = 790.0 ± 190.0,
    GPI_Keq = 0.36 ± 0.11,
    GPI_MW = 63147.0 / 1000,
    PFKP_L = 6.9 ± 4.0,
    PFKP_K_a_F6P = 1.7e-3 ± 0.5e-3,
    PFKP_K_ATP = 250e-6 ± 130e-6,
    PFKP_K_F16BP = 1.5e-3 ± 0.2e-3,
    PFKP_K_ADP = 360e-6 ± 190e-6,
    PFKP_K_Phosphate = 1.1e-3 ± 0.3e-3,
    PFKP_K_a_ADP_reg = 340e-6 ± 100e-6,
    PFKP_K_i_ATP_reg = 1.1e-3 ± 0.3e-3,
    PFKP_K_a_F26BP = 3.4e-7 ± 1.2e-7,
    PFKP_K_i_F26BP = 1.1e-6 ± 0.3e-6,
    PFKP_K_i_Citrate = 3.6e-3 ± 0.3e-3,
    PFKP_Conc = (cell_protein_density / cell_volume_correction) * (8.6e-4 ± 1.2e-4),
    PFKP_Vmax = 80 ± 22,
    PFKP_Keq = 760.0 ± 380.0,
    PFKP_MW = 85596.0 / 1000,
    ALDO_Km_F16BP = 17e-6 ± 5e-6,
    ALDO_Km_DHAP = 430e-6 ± 240e-6,
    ALDO_Kd_DHAP = 12e-6 ± 5e-6,
    ALDO_Km_GAP = 17e-6 ± 14e-6,
    ALDO_Ki_GAP = 1.8e-6 ± 0.7e-6,
    ALDO_Conc = (cell_protein_density / cell_volume_correction) * (3.0e-3 ± 0.3e-3),
    ALDO_Vmax = 29 ± 8,
    ALDO_Keq = 2.7e-6 ± 1.2e-6,
    ALDO_MW = 39420.0 / 1000,
    TPI_Km_DHAP = 470e-6 ± 80e-6,
    TPI_Km_GAP = 11e-6 ± 2e-6,
    TPI_Conc = (cell_protein_density / cell_volume_correction) * (3.0e-3 ± 0.2e-3),
    TPI_Vmax = 6600.0 ± 800.0,
    TPI_Keq = 0.0045 ± 0.0024,
    TPI_MW = 26669.0 / 1000,
    GAPDH_L = 2.9 ± 1.7,
    GAPDH_K_GAP = 1.7e-6 ± 0.1e-6,
    GAPDH_K_a_NAD = 78e-6 ± 9e-6,
    GAPDH_K_i_NAD = 130e-6 ± 34e-6,
    GAPDH_K_a_Phosphate = 1.5e-3 ± 0.4e-3,
    GAPDH_K_i_Phosphate = 5.5e-3 ± 3.8e-3,
    GAPDH_K_BPG = 0.82e-6 ± 0.14e-6,
    GAPDH_K_a_NADH = 12e-6 ± 7e-6,
    GAPDH_K_i_NADH = 2.5e-6 ± 1.1e-6,
    GAPDH_β_i_BPG = 0.22 ± 0.04,
    GAPDH_Conc = (cell_protein_density / cell_volume_correction) * (5.6e-3 ± 0.5e-3),
    GAPDH_Vmax = 130 ± 20,
    # GAPDH_Keq = 15 ± 6,
    GAPDH_Keq = 68 ± 20,
    GAPDH_MW = 36053.0 / 1000,
    PGK_K_BPG = 3e-6 ± 0.7e-6,
    PGK_K_ADP = 42e-6 ± 10e-6,
    PGK_K_ThreePG = 660e-6 ± 220e-6,
    PGK_K_ATP = 580e-6 ± 230e-6,
    PGK_α = 2.1 ± 0.7,
    PGK_β = 0.22 ± 0.11,
    PGK_γ = 2.0 ± 0.8,
    PGK_Conc = (cell_protein_density / cell_volume_correction) * (3.9e-3 ± 0.4e-3),
    PGK_Vmax = 3500.0 ± 900,
    PGK_Keq = 2000 ± 700,
    PGK_MW = 44615.0 / 1000,
    PGAM_Km_ThreePG = 270e-6 ± 60e-6,
    PGAM_Km_TwoPG = 19e-6 ± 7e-6,
    PGAM_Conc = (cell_protein_density / cell_volume_correction) * (2.3e-3 ± 0.2e-3),
    PGAM_Vmax = 3000 ± 900,
    PGAM_Keq = 0.18 ± 0.05,
    PGAM_MW = 28804.0 / 1000,
    ENO_Km_TwoPG = 40e-6 ± 9e-6,
    ENO_Km_PEP = 160e-6 ± 30e-6,
    ENO_Conc = (cell_protein_density / cell_volume_correction) * (9.6e-3 ± 0.5e-3),
    ENO_Vmax = 71 ± 9,
    ENO_Keq = 4.4 ± 1.0,
    ENO_MW = 47169.0 / 1000,
    PKM2_L = 26.0 ± 7.0,
    PKM2_Vmax_a = 540.0 ± 100.0,
    PKM2_Vmax_i = 86.0 ± 41.0,
    PKM2_K_a_PEP = 150e-6 ± 20e-6,
    PKM2_K_ADP = 320e-6 ± 50e-6,
    PKM2_K_Pyruvate = 1800e-6 ± 400e-6,
    PKM2_K_ATP = 940e-6 ± 250e-6,
    PKM2_K_a_F16BP = 1.4e-6 ± 0.4e-6,
    PKM2_K_a_Phenylalanine = 2100e-6 ± 600e-6,
    PKM2_K_i_PEP = 2400e-6 ± 600e-6,
    PKM2_K_i_Phenylalanine = 200e-6 ± 30e-6,
    PKM2_β_i_PEP_ATP = 12.0 ± 5.0,
    PKM2_Keq = 20000.0 ± 6000.0,
    PKM2_Conc = (cell_protein_density / cell_volume_correction) * (7.5e-3 ± 1.0e-3),
    PKM2_MW = 57937.0 / 1000,
    LDH_Km_Pyruvate = 140e-6 ± 20e-6,
    LDH_Km_NADH = 20e-6 ± 2e-6,
    LDH_Km_Lactate = 6.6e-3 ± 0.7e-3,
    LDH_Km_NAD = 220e-6 ± 30e-6,
    LDH_Kd_NADH = 5.7e-6 ± 0.4e-6,
    LDH_Kd_NAD = 810e-6 ± 80e-6,
    LDH_Conc = (cell_protein_density / cell_volume_correction) * (5.1e-3 ± 0.5e-3),
    LDH_Vmax = 320 ± 40,
    LDH_Keq = 13000 ± 5000,
    LDH_MW = 36689.0 / 1000,
    MCT_Km_Lactate = 5.5e-3 ± 0.9e-3,
    MCT_Conc = (cell_protein_density / cell_volume_correction) * (1.5e-4 ± 0.3e-4),
    MCT_Vmax = 470.0 ± 230,
    MCT_Keq = 1.0,
    MCT_MW = 53944.0 / 1000,
    AK_Km_ADP = 0.12e-3,
    AK_Km_ATP = 0.11e-3,
    AK_Km_AMP = 0.09e-3,
    AK_Vmax = 0.03,
    AK_Keq = 0.48,
    AK_MW = 21635.0 / 1000,
    NDPK_Km_ATP = 2e-3,
    NDPK_Km_ADP = 0.1e-3,
    NDPK_Km_NTP = 0.5e-3,
    NDPK_Km_NDP = 0.2e-3,
    NDPK_Vmax = 0.0,
    NDPK_Keq = 2.16, #Average of 2.57 (UDP), 2.71 (GDP) and 1.21 (CTP)
    NDPK_MW = 17149.0 / 1000,
    CK_Km_ATP = 0.16e-3,
    CK_Km_ADP = 0.11e-3,
    CK_Km_Phosphocreatine = 0.9e-3,
    CK_Km_Creatine = 2.3e-3,
    CK_Vmax = 0.0,
    CK_Keq = 0.00598 ± 0.00093,
    CK_MW = 17149.0 / 1000,
    ATPase_Km_ATP = 1e-6,
    ATPase_Km_ADP = 1e-3,
    ATPase_Km_Phosphate = 1e-3,
    ATPase_Keq = 83000.0,
    ATPase_Vmax = 0.0002,
)

glycolysis_params = Measurements.value.(glycolysis_params_w_uncertainty)

#=
    This file contains functions that take M concentrations of glycolytic metabolites as input and
    and generate rates of glycolytic reactions in M/s
=#


"""
    rate(::Enzyme{:GLUT1}, metabs, params)

Calculate rate (M/s units) of GLUT transporter from concentrations (M units) of `Glucose_media` and `Glucose` according to the following equation:

```math
Rate = \\frac{{V_{max} \\cdot Conc}}{{K_{M}^{Glucose}}} \\cdot \\frac{{Glucose_{media} - \\frac{1}{K_{eq}} \\cdot Glucose}}{1 + \\frac{Glucose_{media}}{K_{M}^{Glucose}} + \\frac{Glucose}{K_{M}^{Glucose}}}
```

# Arguments
- `metabs`: LArray or NamedTuple and struct that contains fields Glucose_media, Glucose with corresponding metabolite concentrations of GLUT substrates and products. Glycolysis.jl exports LArray `glycolysis_init_cond` that contains estimates of cellular metabolite concentrations.
- `params`: LArray or NamedTuple and struct of kinetic parameters of GLUT. Glycolysis.jl exports LArray `glycolysis_params` that contains kinetic parameters of all glycolytic enzymes.

# Example
```julia-repl
julia> metabs = (Glucose_media = 25e-3, Glucose = 8e-3,)
julia> Glycolysis.rate_GLUT(metabs, glycolysis_params)
0.02267962645321136
```
"""
@inline function CellMetabolismBase.rate(::Enzyme{:GLUT1, (:Glucose_media,), (:Glucose,)}, metabs, params)
    Rate = (
        (params.GLUT1_Vmax * params.GLUT1_Conc / params.GLUT1_Km_Glucose) *
        (metabs.Glucose_media - (1 / params.GLUT1_Keq) * metabs.Glucose) / (
            1 +
            metabs.Glucose_media / params.GLUT1_Km_Glucose +
            metabs.Glucose / params.GLUT1_Km_Glucose
        )
    )
    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:HK1, (:Glucose, :ATP), (:G6P, :ADP,)}, metabs, params)
    Z = (
        (
            1 +
            (metabs.Glucose / params.HK1_K_Glucose) +
            (metabs.ATP / params.HK1_K_a_ATP) +
            (metabs.G6P / params.HK1_K_G6P) +
            (metabs.G6P / params.HK1_K_a_G6P_cat) +
            (metabs.ADP / params.HK1_K_a_ADP) +
            (params.HK1_β_Glucose_ATP) *
            (metabs.Glucose / params.HK1_K_Glucose) *
            (metabs.ATP / params.HK1_K_a_ATP) +
            (metabs.Glucose / params.HK1_K_Glucose) * (metabs.ADP / params.HK1_K_a_ADP) +
            (metabs.Glucose / params.HK1_K_Glucose) *
            (metabs.G6P / params.HK1_K_a_G6P_cat) +
            (metabs.G6P / params.HK1_K_G6P) * (metabs.ADP / params.HK1_K_a_ADP) +
            (metabs.G6P / params.HK1_K_G6P) * (metabs.G6P / params.HK1_K_a_G6P_cat)
        ) * (1 + (metabs.Phosphate / params.HK1_K_a_Pi)) +
        (1 + (metabs.Glucose / params.HK1_K_Glucose) + (metabs.G6P / params.HK1_K_G6P)) * (metabs.G6P / params.HK1_K_i_G6P_reg)
    )
    Rate =
        (
            (
                params.HK1_Vmax *
                params.HK1_Conc *
                (params.HK1_β_Glucose_ATP) *
                (1 / params.HK1_K_Glucose) *
                (1 / params.HK1_K_a_ATP) *
                (1 + (metabs.Phosphate / params.HK1_K_a_Pi))
            ) * (metabs.Glucose * metabs.ATP - metabs.G6P * metabs.ADP / params.HK1_Keq)
        ) / Z
    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:GPI, (:G6P,), (:F6P,)}, metabs, params)
    Rate = (
        (params.GPI_Vmax * params.GPI_Conc / params.GPI_Km_G6P) *
        (metabs.G6P - (1 / params.GPI_Keq) * metabs.F6P) /
        (1 + metabs.G6P / params.GPI_Km_G6P + metabs.F6P / params.GPI_Km_F6P)
    )
    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:PFKP, (:F6P, :ATP,), (:F16BP, :ADP,)}, metabs, params)

    Z_a_cat = (
        1 +
        (metabs.F6P / params.PFKP_K_a_F6P) +
        (metabs.ATP / params.PFKP_K_ATP) +
        (metabs.F16BP / params.PFKP_K_F16BP) +
        (metabs.ADP / params.PFKP_K_ADP) +
        (metabs.F6P / params.PFKP_K_a_F6P) * (metabs.ATP / params.PFKP_K_ATP) +
        (metabs.F16BP / params.PFKP_K_F16BP) * (metabs.ADP / params.PFKP_K_ADP)
    )
    Z_i_cat = (
        1 +
        (metabs.ATP / params.PFKP_K_ATP) +
        (metabs.F16BP / params.PFKP_K_F16BP) +
        (metabs.ADP / params.PFKP_K_ADP) +
        (metabs.F16BP / params.PFKP_K_F16BP) * (metabs.ADP / params.PFKP_K_ADP)
    )
    Z_a_reg = (
        (1 + metabs.Phosphate / params.PFKP_K_Phosphate) *
        (1 + metabs.ADP / params.PFKP_K_a_ADP_reg) *
        (1 + metabs.F26BP / params.PFKP_K_a_F26BP)
    )
    Z_i_reg = (
        (
            1 +
            metabs.ATP / params.PFKP_K_i_ATP_reg +
            metabs.Phosphate / params.PFKP_K_Phosphate
        ) *
        (1 + metabs.F26BP / params.PFKP_K_i_F26BP) *
        (1 + metabs.Citrate / params.PFKP_K_i_Citrate)
    )

    Rate =
        params.PFKP_Vmax *
        params.PFKP_Conc *
        (metabs.F6P * metabs.ATP - metabs.F16BP * metabs.ADP / params.PFKP_Keq) *
        (1 / params.PFKP_K_a_F6P) *
        (1 / params.PFKP_K_ATP) *
        (Z_a_cat^3) *
        (Z_a_reg^4) /
        ((Z_a_cat^4) * (Z_a_reg^4) + params.PFKP_L * (Z_i_cat^4) * (Z_i_reg^4))

    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:ALDO, (:F16BP,), (:GAP, :DHAP,)}, metabs, params)
    Rate = (
        (params.ALDO_Vmax * params.ALDO_Conc / params.ALDO_Km_F16BP) * (
            (metabs.F16BP - (1 / params.ALDO_Keq) * (metabs.DHAP * metabs.GAP)) / (
                1 +
                metabs.GAP * metabs.DHAP / (params.ALDO_Kd_DHAP * params.ALDO_Km_GAP) +
                metabs.DHAP / params.ALDO_Kd_DHAP +
                metabs.F16BP * metabs.GAP / (params.ALDO_Ki_GAP * params.ALDO_Km_F16BP) +
                metabs.F16BP / params.ALDO_Km_F16BP +
                metabs.GAP * params.ALDO_Km_DHAP /
                (params.ALDO_Kd_DHAP * params.ALDO_Km_GAP)
            )
        )
    )
    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:TPI, (:DHAP,), (:GAP,)}, metabs, params)
    Rate = (
        (params.TPI_Vmax * params.TPI_Conc / params.TPI_Km_DHAP) *
        (metabs.DHAP - (1 / params.TPI_Keq) * metabs.GAP) /
        (1 + (metabs.DHAP / params.TPI_Km_DHAP) + (metabs.GAP / params.TPI_Km_GAP))
    )
    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:GAPDH, (:GAP, :NAD, :Phosphate,), (:BPG, :NADH,)}, metabs, params)
    Z_a =
        (
            1 +
            metabs.GAP / params.GAPDH_K_GAP *
            (1 + metabs.Phosphate / params.GAPDH_K_a_Phosphate) +
            metabs.BPG / params.GAPDH_K_BPG
        ) * (1 + metabs.NAD / params.GAPDH_K_a_NAD + metabs.NADH / params.GAPDH_K_a_NADH)

    Z_i =
        (1 + metabs.NAD / params.GAPDH_K_i_NAD) * (
            1 +
            metabs.GAP / params.GAPDH_K_GAP *
            (1 + metabs.Phosphate / params.GAPDH_K_i_Phosphate) +
            metabs.BPG / params.GAPDH_K_BPG
        ) +
        metabs.NADH / params.GAPDH_K_i_NADH * (
            1 +
            metabs.GAP / params.GAPDH_K_GAP *
            (1 + metabs.Phosphate / params.GAPDH_K_i_Phosphate) +
            params.GAPDH_β_i_BPG * metabs.BPG / params.GAPDH_K_BPG
        )

    Rate = (
        (
            params.GAPDH_Vmax * params.GAPDH_Conc /
            (params.GAPDH_K_GAP * params.GAPDH_K_a_NAD * params.GAPDH_K_a_Phosphate)
        ) *
        Z_a^3 *
        (
            metabs.GAP * metabs.NAD * metabs.Phosphate -
            (1 / params.GAPDH_Keq) * metabs.BPG * metabs.NADH
        ) / (Z_a^4 + params.GAPDH_L * Z_i^4)
    )
    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:PGK, (:BPG, :ADP,), (:ThreePG, :ATP,)}, metabs, params)
    Rate = (
        (
            params.PGK_Vmax * params.PGK_Conc /
            (params.PGK_α * params.PGK_K_BPG * params.PGK_K_ADP)
        ) * (
            metabs.BPG * metabs.ADP - (1 / params.PGK_Keq) * (metabs.ThreePG * metabs.ATP)
        ) / (
            1 +
            metabs.BPG / params.PGK_K_BPG +
            metabs.ADP / params.PGK_K_ADP +
            metabs.ThreePG / params.PGK_K_ThreePG +
            metabs.ATP / params.PGK_K_ATP +
            metabs.BPG * metabs.ADP / (params.PGK_α * params.PGK_K_BPG * params.PGK_K_ADP) +
            metabs.ThreePG * metabs.ATP /
            (params.PGK_β * params.PGK_K_ThreePG * params.PGK_K_ATP) +
            metabs.ThreePG * metabs.ADP /
            (params.PGK_γ * params.PGK_K_ThreePG * params.PGK_K_ADP)
        )
    )
    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:PGAM, (:ThreePG,), (:TwoPG,)}, metabs, params)
    Rate = (
        (params.PGAM_Vmax * params.PGAM_Conc / params.PGAM_Km_ThreePG) *
        (metabs.ThreePG - (1 / params.PGAM_Keq) * metabs.TwoPG) / (
            1 +
            metabs.ThreePG / params.PGAM_Km_ThreePG +
            metabs.TwoPG / params.PGAM_Km_TwoPG
        )
    )
    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:ENO, (:TwoPG,), (:PEP,)}, metabs, params)
    Rate = (
        (params.ENO_Vmax * params.ENO_Conc / params.ENO_Km_TwoPG) *
        (metabs.TwoPG - (1 / params.ENO_Keq) * metabs.PEP) /
        (1 + metabs.TwoPG / params.ENO_Km_TwoPG + metabs.PEP / params.ENO_Km_PEP)
    )
    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:PKM2, (:PEP, :ADP,), (:Pyruvate, :ATP,)}, metabs, params)

    Z_a_cat = (
        1 +
        (metabs.PEP / params.PKM2_K_a_PEP) +
        (metabs.ATP / params.PKM2_K_ATP) +
        (metabs.ADP / params.PKM2_K_ADP) +
        (metabs.PEP / params.PKM2_K_a_PEP) * (metabs.ADP / params.PKM2_K_ADP) +
        (metabs.Pyruvate / params.PKM2_K_Pyruvate) * (metabs.ATP / params.PKM2_K_ATP) +
        (metabs.PEP / params.PKM2_K_a_PEP) * (metabs.ATP / params.PKM2_K_ATP) +
        (metabs.Pyruvate / params.PKM2_K_Pyruvate) * (metabs.ADP / params.PKM2_K_ADP)
    )
    Z_i_cat = (
        1 +
        (metabs.PEP / params.PKM2_K_i_PEP) +
        (metabs.ATP / params.PKM2_K_ATP) +
        (metabs.ADP / params.PKM2_K_ADP) +
        (metabs.PEP / params.PKM2_K_i_PEP) * (metabs.ADP / params.PKM2_K_ADP) +
        (metabs.Pyruvate / params.PKM2_K_Pyruvate) * (metabs.ATP / params.PKM2_K_ATP) +
        params.PKM2_β_i_PEP_ATP *
        (metabs.PEP / params.PKM2_K_i_PEP) *
        (metabs.ATP / params.PKM2_K_ATP) +
        (metabs.Pyruvate / params.PKM2_K_Pyruvate) * (metabs.ADP / params.PKM2_K_ADP)
    )
    Z_a_reg = (
        (1 + metabs.F16BP / params.PKM2_K_a_F16BP) *
        (1 + metabs.Phenylalanine / params.PKM2_K_a_Phenylalanine)
    )
    Z_i_reg = (1 + metabs.Phenylalanine / params.PKM2_K_i_Phenylalanine)

    Rate =
        params.PKM2_Conc *
        (metabs.ADP * metabs.PEP - metabs.ATP * metabs.Pyruvate / params.PKM2_Keq) *
        (
            (params.PKM2_Vmax_a * (1.0 / params.PKM2_K_a_PEP) * (1.0 / params.PKM2_K_ADP)) *
            (Z_a_cat^3) *
            (Z_a_reg^4) +
            params.PKM2_L *
            (params.PKM2_Vmax_i * (1.0 / params.PKM2_K_i_PEP) * (1.0 / params.PKM2_K_ADP)) *
            (Z_i_cat^3) *
            (Z_i_reg^4)
        ) / ((Z_a_cat^4) * (Z_a_reg^4) + params.PKM2_L * (Z_i_cat^4) * (Z_i_reg^4))

    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:LDH, (:Pyruvate, :NADH,), (:Lactate, :NAD,)}, metabs, params)
    Rate = (
        (
            params.LDH_Vmax * params.LDH_Conc /
            (params.LDH_Km_Pyruvate * params.LDH_Kd_NADH)
        ) * (
            metabs.Pyruvate * metabs.NADH -
            (1 / params.LDH_Keq) * (metabs.Lactate * metabs.NAD)
        ) / (
            1 +
            metabs.Pyruvate * params.LDH_Km_NADH /
            (params.LDH_Kd_NADH * params.LDH_Km_Pyruvate) +
            metabs.Lactate * params.LDH_Km_NAD /
            (params.LDH_Kd_NAD * params.LDH_Km_Lactate) +
            metabs.NADH / params.LDH_Kd_NADH +
            metabs.Lactate * metabs.NAD / (params.LDH_Kd_NAD * params.LDH_Km_Lactate) +
            metabs.Lactate * metabs.NADH * params.LDH_Km_NAD /
            (params.LDH_Kd_NAD * params.LDH_Kd_NADH * params.LDH_Km_Lactate) +
            metabs.Pyruvate * metabs.NADH / (params.LDH_Kd_NADH * params.LDH_Km_Pyruvate) +
            metabs.NAD / params.LDH_Kd_NAD +
            metabs.Pyruvate * metabs.NAD * params.LDH_Km_NADH /
            (params.LDH_Kd_NAD * params.LDH_Kd_NADH * params.LDH_Km_Pyruvate)
        )
    )
    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:MCT, (:Lactate,), (:Lactate_media,)}, metabs, params)
    Rate = (
        (params.MCT_Vmax * params.MCT_Conc / params.MCT_Km_Lactate) *
        (metabs.Lactate - (1 / params.MCT_Keq) * metabs.Lactate_media) / (
            1 +
            metabs.Lactate / params.MCT_Km_Lactate +
            metabs.Lactate_media / params.MCT_Km_Lactate
        )
    )
    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:AK, (:ADP, :ADP,), (:ATP, :AMP,)}, metabs, params)
    Rate = (
        (params.AK_Vmax / (params.AK_Km_ADP^2)) *
        (metabs.ADP^2 - (1 / params.AK_Keq) * (metabs.ATP * metabs.AMP)) / (
            (1 + metabs.ADP / params.AK_Km_ADP + metabs.ATP / params.AK_Km_ATP) *
            (1 + metabs.ADP / params.AK_Km_ADP + metabs.AMP / params.AK_Km_AMP)
        )
    )
    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:NDPK, (:ATP, :NDP,), (:ADP, :NTP,)}, metabs, params)
    Rate = (
        (params.NDPK_Vmax / (params.NDPK_Km_ATP * params.NDPK_Km_NDP)) * (
            (metabs.ATP * metabs.NDP - (1 / params.NDPK_Keq) * (metabs.NTP * metabs.ADP)) / (
                (1 + metabs.ATP / params.NDPK_Km_ATP + metabs.ADP / params.NDPK_Km_ADP) *
                (1 + metabs.NTP / params.NDPK_Km_NTP + metabs.NDP / params.NDPK_Km_NDP)
            )
        )
    )
    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:CK, (:ATP, :Creatine,), (:Phosphocreatine, :ADP,)}, metabs, params)
    Rate = (
        (params.CK_Vmax / (params.CK_Km_ATP * params.CK_Km_Creatine)) * (
            (
                metabs.ATP * metabs.Creatine -
                (1 / params.CK_Keq) * (metabs.Phosphocreatine * metabs.ADP)
            ) / (
                (1 + metabs.ATP / params.CK_Km_ATP + metabs.ADP / params.CK_Km_ADP) * (
                    1 +
                    metabs.Phosphocreatine / params.CK_Km_Phosphocreatine +
                    metabs.Creatine / params.CK_Km_Creatine
                )
            )
        )
    )
    return Rate
end


@inline function CellMetabolismBase.rate(::Enzyme{:ATPase, (:ATP,), (:ADP, :Phosphate,)}, metabs, params)
    Rate =
        (params.ATPase_Vmax / params.ATPase_Km_ATP) *
        (
            1 / (
                1 +
                metabs.ATP / params.ATPase_Km_ATP +
                metabs.ADP / params.ATPase_Km_ADP +
                metabs.Phosphate / params.ATPase_Km_Phosphate +
                (metabs.ADP / params.ATPase_Km_ADP) *
                (metabs.Phosphate / params.ATPase_Km_Phosphate)
            )
        ) *
        (metabs.ATP - metabs.ADP * metabs.Phosphate / params.ATPase_Keq)
    return Rate
end
