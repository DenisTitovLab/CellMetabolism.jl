using CellMetabolism, CellMetabolismBase, OrdinaryDiffEqFIRK
using CairoMakie, Dates

# Constants
const F26BP_RANGE = (2e-9, 0.5e-3)
const F26BP_PHYS_RANGE = (1e-8, 1e-4)

# Total ATP supply capacity (respiration + glycolysis) - kept constant
# 0.3 is because glycolysis model cannot go faster than ~30% of HK1_Vmax * HK1_Conc
# 2 is because glycolysis produces 2 ATP per HK1 reaction flux
const TOTAL_ATP_SUPPLY_VMAX =
    0.3 * 2 *
    CellMetabolism.glycolysis_params.HK1_Vmax *
    CellMetabolism.glycolysis_params.HK1_Conc

# ATPase - kept constant
const ATPASE_VMAX = 0.2 * TOTAL_ATP_SUPPLY_VMAX

base_params_nt = NamedTuple(CellMetabolism.glycolysis_params)
f26bp_values = logrange(F26BP_RANGE..., length = 100)

# List of all glycolytic enzyme concentration parameters to scale
# ATPase is excluded as it has no concentration parameter
const GLYCOLYTIC_ENZYME_CONC_PARAMS = (
    :GLUT1_Conc,
    :HK1_Conc,
    :GPI_Conc,
    :PFKP_Conc,
    :ALDO_Conc,
    :TPI_Conc,
    :GAPDH_Conc,
    :PGK_Conc,
    :PGAM_Conc,
    :ENO_Conc,
    :PKM2_Conc,
    :LDH_Conc,
    :MCT_Conc,
    :AK_Conc,
    :NDPK_Conc,
)

"""
    scale_glycolysis_enzyme_concentrations(base_params, multiplier)

Scale all glycolytic enzyme concentrations by the given multiplier.
This properly adjusts glycolysis Vmax by scaling enzyme concentrations
rather than just modifying a single enzyme's Vmax.

ATPase is excluded as it has no concentration parameter.
"""
function scale_glycolysis_enzyme_concentrations(base_params::NamedTuple, multiplier::Real)
    scaled_params = Dict{Symbol,Any}()
    for (key, val) in pairs(base_params)
        if key in GLYCOLYTIC_ENZYME_CONC_PARAMS
            scaled_params[key] = val * multiplier
        else
            scaled_params[key] = val
        end
    end
    return NamedTuple{Tuple(keys(scaled_params))}(values(scaled_params))
end

# Setup the pathway - use proper enzyme definitions with regulators matching CellMetabolism
modified_glycolysis_pathway = MetabolicPathway(
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
        (:ATPase, (:ATP,), (:ADP, :Phosphate)),
        (:ATPSynthase, (:ADP, :Phosphate), (:ATP,)),
    ),
)

@inline function CellMetabolismBase.rate(
    ::Enzyme{:ATPSynthase,(:ADP, :Phosphate),(:ATP,)},
    metabs,
    params,
)
    ATPSynthase_Km_ATP = 1000
    ATPSynthase_Km_ADP = 1e-5
    ATPSynthase_Km_Phosphate = 1e-5

    Rate =
        (params.ATPSynthase_Vmax / (ATPSynthase_Km_ADP * ATPSynthase_Km_Phosphate)) *
        (
            1 / (
                1 +
                metabs.ATP / ATPSynthase_Km_ATP +
                metabs.ADP / ATPSynthase_Km_ADP +
                metabs.Phosphate / ATPSynthase_Km_Phosphate +
                (metabs.ADP / ATPSynthase_Km_ADP) *
                (metabs.Phosphate / ATPSynthase_Km_Phosphate)
            )
        ) *
        (metabs.ADP * metabs.Phosphate - metabs.ATP / params.ATPSynthase_Keq)
    return Rate
end

"""
    run_simulation(glycolysis_multiplier, respiration_vmax, atpase_vmax, f26bp_values)

Run simulation with given respiration Vmax and glycolysis multiplier.

The glycolysis_multiplier scales all glycolytic enzyme concentrations.
Since glycolysis produces 2 ATP per glucose, baseline glycolysis ATP production
rate = 2 * HK1_Conc * HK1_Vmax. After scaling by multiplier:
    glycolysis_vmax = 2 * (HK1_Conc * multiplier) * HK1_Vmax

To maintain respiration_vmax + glycolysis_vmax = TOTAL_ATP_SUPPLY_VMAX:
    glycolysis_multiplier = (1 - respiration_fraction)
"""
function run_simulation(glycolysis_multiplier, respiration_vmax, atpase_vmax, f26bp_values)
    # Scale all enzyme concentrations by the glycolysis multiplier
    scaled_params_nt = scale_glycolysis_enzyme_concentrations(base_params_nt, glycolysis_multiplier)

    params = CellMetabolism.LabelledArrays.LVector(
        merge(
            scaled_params_nt,
            (
                ATPSynthase_Vmax = respiration_vmax,
                ATPSynthase_Keq = 1_000_000.0,
                ATPase_Vmax = atpase_vmax,
            ),
        ),
    )

    # Generate different initial conditions for each F26BP value
    init_conditions = Vector{typeof(CellMetabolism.glycolysis_init_conc)}()
    for f26bp in f26bp_values
        init_cond = deepcopy(CellMetabolism.glycolysis_init_conc)
        init_cond.F26BP = f26bp
        push!(init_conditions, init_cond)
    end

    # Run to steady state
    tspan = (0.0, 1e8)
    ensemble_prob = CellMetabolismBase.make_EnsembleProblem(
        modified_glycolysis_pathway,
        init_conditions,
        params;
        tspan = tspan,
    )
    ensemble_sol = solve(
        ensemble_prob,
        RadauIIA9(),
        trajectories = length(init_conditions),
        abstol = 1e-15,
        reltol = 1e-8,
    )

    # Calculate enzyme rates
    enzyme_rate_vals = []
    for (i, f26bp) in enumerate(f26bp_values)
        push!(
            enzyme_rate_vals,
            [
                CellMetabolismBase.rates(
                    modified_glycolysis_pathway,
                    ensemble_sol[i].u[end],
                    params,
                )...,
            ],
        )
    end
    enzyme_rate_matrix = hcat(enzyme_rate_vals...)
    enzyme_names = [CellMetabolismBase.enzymes(modified_glycolysis_pathway)...]
    enzyme_rates_dict =
        Dict(enzyme_names[i] => enzyme_rate_matrix[i, :] for i in eachindex(enzyme_names))

    return enzyme_rates_dict
end

# Define conditions with varying respiration/glycolysis ratio over 10-fold range
# respiration_fraction: fraction of total ATP supply from respiration
# glycolysis uses the remainder
#
# Key constraint: respiration_vmax + glycolysis_vmax = TOTAL_ATP_SUPPLY_VMAX
# Since glycolysis produces 2 ATP per glucose, glycolysis_multiplier = (1 - frac)
respiration_fractions = [0.2, 0.5, 0.8]  # 10-fold range: 0.1/0.9 to 0.9/0.1
conditions = [
    (
        resp_vmax = TOTAL_ATP_SUPPLY_VMAX * frac,
        glyc_multiplier = 1 - frac,
        title = let r = round(frac*10, digits=1), g = round((1-frac)*10, digits=1)
            "Resp:Glyc = $(isinteger(r) ? Int(r) : r):$(isinteger(g) ? Int(g) : g)"
        end,
        col = i,
    )
    for (i, frac) in enumerate(respiration_fractions)
]

# Run simulations for all conditions
all_results = Dict{Int,Dict}()
for cond in conditions
    all_results[cond.col] = run_simulation(cond.glyc_multiplier, cond.resp_vmax, ATPASE_VMAX, f26bp_values)
end

# Plotting configuration - unified enzyme definitions
enzymes = [
    (sym = :PKM2, label = "Glycolysis", color = Makie.RGB(0xFF/255, 0x2F/255, 0x92/255), lw = 6),
    (sym = :ATPSynthase, label = "Respiration", color = Makie.RGB(0, 0, 0), lw = 6),
    (sym = :ATPase, label = "ATPase", color = Makie.RGB(0x00/255, 0x91/255, 0x93/255), lw = 3),
]

# Compute common y-axis limits across all conditions
ymax_raw = maximum(maximum(all_results[c.col][e.sym]) for c in conditions for e in enzymes)
ymin_raw = minimum(minimum(all_results[c.col][e.sym]) for c in conditions for e in enzymes)
ymax = ymax_raw + 0.05 * ymax_raw
ymin = 0 - 0.1 * ymax_raw

# Create the side-by-side plot
set_theme!()
inch = 96
fig = Figure(size = (6.5inch, 3inch))

# Loop over conditions to create each panel
first_ax = nothing
for cond in conditions
    rates = all_results[cond.col]

    is_first = cond.col == 1
    ax = Axis(
        fig[1, cond.col],
        xlabel = "[F26BP], M",
        ylabel = "Rate of ATP cons. or prod., M/s",
        title = cond.title,
        xscale = log10,
        xticklabelsize = 10,
        yticklabelsize = 10,
        yticklabelrotation = π/2,
        limits = (F26BP_RANGE, (ymin, ymax)),
    )

    # Link y-axis and hide decorations for non-first axes
    if is_first
        global first_ax = ax
    else
        linkyaxes!(first_ax, ax)
        hideydecorations!(ax, grid = false)
    end

    # Add physiological F26BP shading
    vspan!(ax, F26BP_PHYS_RANGE..., color = (:lightblue, 0.3))
    text!(
        ax,
        sqrt(prod(F26BP_PHYS_RANGE)),
        0 - 0.075 * (ymax - ymin),
        text = "Physiological [F26BP]",
        color = :dodgerblue,
        align = (:center, :bottom),
        fontsize = 10,
    )

    # Plot enzymes
    for enz in enzymes
        lines!(ax, f26bp_values, rates[enz.sym], color = enz.color, linewidth = enz.lw)
    end

    # Add inline labels with custom positions
    label_offset = 0.05 * (ymax - ymin)
    left_idx = 5  # index near left side
    right_idx = length(f26bp_values) - 5  # index near right side

    # Glycolysis: left side, slightly below line
    text!(ax, f26bp_values[left_idx], rates[:PKM2][left_idx] - label_offset,
          text = "Glycolysis", color = enzymes[1].color,
          align = (:left, :center), fontsize = 12)

    # Respiration: left side, slightly above line
    text!(ax, f26bp_values[left_idx], rates[:ATPSynthase][left_idx] + label_offset,
          text = "Respiration", color = enzymes[2].color,
          align = (:left, :center), fontsize = 12)

    # ATPase: right side, slightly below line
    text!(ax, f26bp_values[right_idx], rates[:ATPase][right_idx] - label_offset,
          text = "ATPase", color = enzymes[3].color,
          align = (:right, :center), fontsize = 12)
end

# Add A, B, C labels to top right corner of each plot
for (i, label) in enumerate(["A", "B", "C"])
    Label(fig[1, i, TopLeft()], label, fontsize = 12, padding = (0, 0, 5, 0), halign = :right)
end

display(fig)

# save("Plots/$(Dates.format(today(), "yyyy_mm_dd"))_Effect_of_F26BP_varying_resp_glyc_ratio.pdf", fig)
