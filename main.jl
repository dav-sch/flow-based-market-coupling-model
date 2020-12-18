using Pkg
using CSV, XLSX, DataFrames
using Gurobi
using JuMP
using Statistics
using Plots
using DataFrames

## LOAD DATA, FUNCTIONS, COMPUTE PRELIMINARIES
#cd("...") #may have to change
include("data_and_functions.jl")
include("fbmc_preliminaries.jl")

## PARAMETER CHOICES
# Flow-based preliminaries
cne_alpha = 0.1
gsk_cne = gsk_flat_unit
gsk_mc = gsk_flat_unit
frm = 0.2
include_cb_lines = true
include("cne_selection.jl")

# Market coupling
max_ntc = 1000
cost_curt_mc = 0
np_within_fbmc_0 = false
exports_abroad_0 = true

## MODELING
include("model_D-2_base_case.jl")
include("model_D-1_market_coupling.jl")
include("model_D-0_congestion_management.jl")

## EVALUATION
include("eval_results.jl")
