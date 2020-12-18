### Set path:
cd("./data")

### Import data:
df_bus_load = CSV.read("df_bus_load_added_abroad_final.csv")
df_bus = CSV.read("df_bus_final.csv")
df_branch = CSV.read("df_branch_final.csv")
df_plants = CSV.read("df_gen_final.csv", copycols=true)
#df_plants = CSV.read("df_gen_final_high_RES.csv", copycols=true)
incidence = CSV.read("matrix_A_final.csv")
susceptance = CSV.read("matrix_Bd_final.csv")

xf_renew = XLSX.readxlsx("data_renew_2015.xlsx")
df_pv = DataFrame(xf_renew["pv"][:][2:end,:])
rename!(df_pv, Dict(names(df_pv)[i] => Symbol.(xf_renew["pv"][:][1,:])[i] for i = 1:ncol(df_pv)))
df_wind = DataFrame(xf_renew["onshore"][:][2:end,:])
rename!(df_wind, Dict(names(df_wind)[i] => Symbol.(xf_renew["onshore"][:][1,:])[i] for i = 1:ncol(df_wind)))
df_wind_off = DataFrame(xf_renew["offshore"][:][2:end,:])
rename!(df_wind_off, Dict(names(df_wind_off)[i] => Symbol.(xf_renew["offshore"][:][1,:])[i] for i = 1:ncol(df_wind_off)))

# Adjustment of capacities:
df_branch.Pmax = 0.75*df_branch.Pmax
df_bus.ZoneRes = df_bus.Zone


### Sets:
## General sets:
T = 1:size(df_bus_load,1)
R = ["PV","Wind","Wind Offshore"]
P = df_plants.GenID[[!(x in R) for x in df_plants[:,:Type]]]
N = df_bus[:,:BusID]
L = df_branch[:,:BranchID]
Z = sort(unique(df_bus[:,:Zone]))

# This function is used to assign units to zones,
# in case new zone configurations are used in df_bus
function replaced_zones()
	zone_p_new = []
	for i in df_plants.OnBus
		zone_p_new = vcat(zone_p_new,df_bus.Zone[df_bus[:,:BusID].==i][1])
		#df_plants.Zone[df_plants[:,:GenID].==p] = zone_p_new
	end
	return zone_p_new
end

df_plants.Zone = replaced_zones()

## Flow-based sets:
Z_FBMC = Z[1:(length(Z)-3)]
Z_not_in_FBMC = Z[(length(Z)-2):length(Z)]
N_FBMC = df_bus.BusID[[x in Z_FBMC for x in df_bus[:,:Zone]]]
N_not_in_FBMC = df_bus.BusID[[!(x in Z_FBMC) for x in df_bus[:,:Zone]]]

## Redispatch:
#P_RD = df_plants.GenID[[x in ["Hard Coal", "Gas/CCGT"] for x in df_plants[:,:Type]] .&
#					   [x in ["1","2","3","Import/Export_1","Import/Export_2","Import/Export_3"] for x in df_plants[:,:Zone]]]
P_RD = df_plants.GenID[[x in ["Hard Coal", "Gas/CCGT"] for x in df_plants[:,:Type]] .&
					   [x in Z_FBMC for x in df_plants[:,:Zone]]]

## Mapping:
n_in_z = Dict(map(z -> z => [n for n in N if df_bus[df_bus[:,:BusID].==n, :Zone][1] == z], Z))
p_at_n = Dict(map(n -> n => [p for p in P if df_plants[df_plants[:,:GenID].==p, :OnBus][1] == n], N))
p_rd_at_n = Dict(map(n-> n=> [p for p in P_RD if df_plants[df_plants[:,:GenID].==p, :OnBus][1] == n], N))
p_in_z = Dict(map(z -> z => [p for p in P if df_plants[df_plants[:,:GenID].==p, :Zone][1] == z], Z))

z_to_z = Dict(map(z-> z=> [zz for zz in Z if
	     (zz in df_bus.Zone[[(x in df_branch.FromBus[(x in n_in_z[z] for x in df_branch.FromBus) .|
					              (x in n_in_z[z] for x in df_branch.ToBus)]) .|
	 		 (x in df_branch.ToBus[(x in n_in_z[z] for x in df_branch.FromBus) .|
				 	              (x in n_in_z[z] for x in df_branch.ToBus)]) for x in df_bus.BusID]]) .& (zz!=z)], Z))

### Calculations and functions:
## Susceptance matrices:
line_sus_mat = convert(Matrix, susceptance)*convert(Matrix, incidence)
node_sus_mat = transpose(convert(Matrix, incidence))*convert(Matrix, susceptance)*convert(Matrix, incidence)

function get_line_sus(l,n)
	return line_sus_mat[findfirst(L .== l), findfirst(N .== n)]
end

function get_node_sus(n,nn)
	return node_sus_mat[findfirst(N .== n), findfirst(N .== nn)]
end

H_mat = Dict((l,n) => get_line_sus(l,n)
	for (l,l) in enumerate(L), (n,n) in enumerate(N))

B_mat = Dict((n,nn) => get_node_sus(n,nn)
	for (n,n) in enumerate(N), (nn,nn) in enumerate(N))

## Marginal costs:
function get_mc(p)
	return df_plants[df_plants[:,:GenID].==p, :Costs][1]
end

function find_maximum_mc()
	max_temp = 0
    for p in P_RD
        mc_temp = get_mc(p)
        if mc_temp > max_temp
            max_temp = mc_temp
        end
    end
	return max_temp
end

## Demand:
function get_dem(t,n)
	return df_bus_load[t,Symbol.(n)]
end

## Renewables:
function create_res_table()
	res_temp = zeros(Float64, length(T), length(N), length(R))
	for n in N, r in R
		zone_temp = df_bus.ZoneRes[df_bus[:,:BusID].==n][1]
		cap_temp = sum(df_plants.Pmax[(df_plants[:,:Type].==r) .&
		                              (df_plants[:,:OnBus].==n)])
		if r == "PV"
			share_temp = df_pv[:,Symbol.(zone_temp)]
		elseif r == "Wind"
			share_temp = df_wind[:,Symbol.(zone_temp)]
		else
			share_temp = df_wind_off[:,Symbol.(zone_temp)]
		end
		res_temp[:, findfirst(N .== n), findfirst(R .== r)] = cap_temp*share_temp
	end
	return res_temp
end

res_table = create_res_table()

function get_renew(t,n)
	return sum(res_table[findfirst(T .== t), findfirst(N .== n), findfirst(R .== r)] for r in R)
end

## Get conventional capacity
function get_gen_up(p)
	return df_plants.Pmax[df_plants[:,:GenID].==p][1]
end

## Get line capacity
function get_line_cap(l)
	return df_branch.Pmax[df_branch[:,:BranchID].==l][1]
end

## Get line capacity
function find_cross_border_lines()
	cb_lines_temp = []
	for l in L
		from_zone_temp = df_bus.Zone[df_bus[:,:BusID].==df_branch.FromBus[df_branch[:,:BranchID].==l][1]][1]
		to_zone_temp = df_bus.Zone[df_bus[:,:BusID].==df_branch.ToBus[df_branch[:,:BranchID].==l][1]][1]
		if (from_zone_temp in Z_FBMC && to_zone_temp in Z_FBMC && (from_zone_temp!=to_zone_temp))
			cb_lines_temp = vcat(cb_lines_temp,l)
		end
	end
	return cb_lines_temp
end

cd("..")
