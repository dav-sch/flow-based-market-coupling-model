cd("./results/D-1_market_coupling")
d_1_gen_costs = convert(Array{Float64,2}, CSV.read(string("df_d_1_gen_costs.csv")))
d_1_nodal_price= convert(Array{Float64,2}, CSV.read(string("df_d_1_nodal_price.csv")))
cd("..")
cd("..")

cd("./results/D-0_congestion_management")
d_0_curt = convert(Array{Float64,2}, CSV.read(string("df_d_0_curt.csv")))
d_0_rd_pos = convert(Array{Float64,2}, CSV.read(string("df_d_0_rd_pos.csv")))
d_0_rd_neg = convert(Array{Float64,2}, CSV.read(string("df_d_0_rd_neg.csv")))
d_0_redispatch_costs = convert(Array{Float64,2}, CSV.read(string("df_d_0_redispatch_costs.csv")))
cd("..")
cd("..")

function count_unique(v::Array{Float64,1})
    s = length(unique(v))
    return s
end

function price_zones(a::Array{Float64,2})
    v_temp = vec(convert(Array, mapcols!(count_unique, DataFrame(transpose(a[:,1:3])))))
    df = DataFrame(OnePriceZone = sum(v_temp .== 1),
                   TwoPriceZones = sum(v_temp .== 2),
                   ThreePriceZones = sum(v_temp .== 3))
    return df
end

## D-1 MARKET COUPLING
# Generation costs:
d_1_costs_temp = round.(reshape(sum(d_1_gen_costs[:,1:3],dims=1),:)/1e6, digits=3)

# Price zones:
price_zones_temp = vec(convert(Array,price_zones(d_1_nodal_price)))

## D-0 CONGESTION MANAGEMENT
# Amounts:
d_0_rd_pos_temp = round(sum(d_0_rd_pos)/1e3, digits=3)
d_0_rd_neg_temp = round(sum(d_0_rd_neg)/1e3, digits=3)
d_0_curt_temp = round(sum(d_0_curt)/1e3, digits=3)
#round(sum(d_1_curt)/1e3, digits=3)

# Costs:
d_0_redispatch_costs_temp = round(sum(d_0_redispatch_costs)/1e6, digits=3)
costs_curtailment_eval = 100
d_0_curt_costs_temp = round(sum(d_0_curt)*costs_curtailment_eval/1e6, digits=3)

## PRINT:
println("-------------------------------------------------------------------")
println(" \t \t Zone 1 \t Zone 2 \t Zone 3")
println("D-1 Gen. costs \t ", d_1_costs_temp[1], "\t ", d_1_costs_temp[2], "\t ", d_1_costs_temp[3])
println("-------------------------------------------------------------------")
println(" \t \t 1 zone \t 2 zones \t 3 zones")
println("Price zones [h]  ", price_zones_temp[1], "\t \t ", price_zones_temp[2], "\t \t ", price_zones_temp[3])
println("-------------------------------------------------------------------")
println(" \t \t Pos. RD \t Neg. RD \t Curt.")
println("Amounts [GWh] \t ", d_0_rd_pos_temp, "\t \t ", d_0_rd_neg_temp, "\t \t ", d_0_curt_temp)
println("-------------------------------------------------------------------")
println(" \t \t Redispatch \t Curtailment")
println("Costs [mi. EUR]  ", d_0_redispatch_costs_temp, "\t \t ", d_0_curt_costs_temp)
println("-------------------------------------------------------------------")
