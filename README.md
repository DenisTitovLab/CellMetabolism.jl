# CellMetabolism

[![Build Status](https://github.com/DenisTitovLab/CellMetabolism.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DenisTitovLab/CellMetabolism.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/DenisTitovLab/CellMetabolism.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DenisTitovLab/CellMetabolism.jl)
[![JET](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Overview

CellMetabolism.jl is a package for simulation and analysis of human cell metabolism using a experimentally determined enzyme rate equations. The goal of the package is to provide a convenient interface to automatically convert a list of enzymes into a form that can be used by ordinary differential equations (ODE) solvers of [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) and [Scientific Machine Learning (SciML)](https://sciml.ai) ecosystems in [Julia Programming Language](https://julialang.org). CellMetabolism.jl is powered by [CellMetabolismBase.jl](https://github.com/DenisTitovLab/CellMetabolismBase.jl) that enables simulation and analysis of any metabolic pathway using user-provided enzyme rate equations.  

Currently, CellMetabolism.jl supports the glycolysis pathway. The longterm goal is to support all of the human energy metabolism and biosynthesis pathways to allow simulation of intracellular metabolism of any human cell type. Examples of the analyses that can be performed can be found in the [examples](https://github.com/DenisTitovLab/CellMetabolism.jl/tree/main/examples) folder and our [recent publication](https://doi.org/10.1016/j.bpj.2025.03.037).

## Features

- Simulate glycolysis activity using [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) at a wide range of parameter values and initial conditions
- Estimate model prediction uncertainty for comparison with experimental data
- Perform validation to ensure that enzyme rate equations and ODEs are assembled correctly

## Roadmap

- Add support for additional metabolic pathways, including TCA cycle, pentose phosphate pathway, and fatty acid metabolism
- Add additional functionality through updates of [CellMetabolismBase.jl](https://github.com/DenisTitovLab/CellMetabolismBase.jl) such as support for Global Sensitivity Analysis (GSA) to identify parameters that control specific pathway behaviour, ability to simulate isotope tracing equations, support for units of parameters and initial conditions to ensure that the models are dimensionally consistent, etc.