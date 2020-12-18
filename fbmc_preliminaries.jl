MWBase = 380^2
slack_node = 68
slack_position = findfirst(N .== slack_node)

# Build nodal PTDFs
line_sus_mat = convert(Matrix, susceptance)/MWBase*convert(Matrix, incidence)
node_sus_mat = transpose(convert(Matrix, incidence))*convert(Matrix, susceptance)/MWBase*convert(Matrix, incidence)

line_sus_mat_ = line_sus_mat[:, 1:end .!= slack_position]
node_sus_mat_ = node_sus_mat[1:end .!= slack_position, 1:end .!= slack_position]

PTDF = line_sus_mat_*inv(node_sus_mat_)
zero_column = zeros(Float64, length(L), 1)
PTDF = hcat(PTDF[:,1:(slack_position-1)], zero_column, PTDF[:,slack_position:end])

PTDF = PTDF[:,findall(x->x in N_FBMC, N)]

println("Nodal PTDF built.")

# Build flat GSK
function get_gsk_flat()
	gsk_temp = zeros(Float64, length(N_FBMC), length(Z_FBMC))
	for n in N_FBMC
		zone_temp = df_bus.Zone[df_bus[:,:BusID].==n][1]
		gsk_value_temp = 1/size(df_bus[df_bus[:,:Zone].==zone_temp,:],1)
		gsk_temp[findfirst(N_FBMC .== n), findfirst(Z_FBMC .== zone_temp)] = gsk_value_temp
	end
	return gsk_temp
end

gsk_flat = get_gsk_flat()
sum(gsk_flat, dims=1)

function get_gsk_flat_unit()
	gsk_temp = zeros(Float64, length(N_FBMC), length(Z_FBMC))
	for n in N_FBMC
		zone_temp = df_bus.Zone[df_bus[:,:BusID].==n][1]
		conv_nodes_in_zone = unique(df_plants.OnBus[[x == zone_temp for x in df_plants[:,:Zone]] .& [x in P for x in df_plants[:,:GenID]]])
		if n in conv_nodes_in_zone
			gsk_value_temp = 1/length(conv_nodes_in_zone)
			gsk_temp[findfirst(N_FBMC .== n), findfirst(Z_FBMC .== zone_temp)] = gsk_value_temp
		end
	end
	return gsk_temp
end

gsk_flat_unit = get_gsk_flat_unit()
sum(gsk_flat_unit, dims=1)

function get_gsk_pmax()
	gsk_temp = zeros(Float64, length(N_FBMC), length(Z_FBMC))
	for n in N_FBMC
		zone_temp = df_bus.Zone[df_bus[:,:BusID].==n][1]
		conv_nodes_in_zone = unique(df_plants.OnBus[[x == zone_temp for x in df_plants[:,:Zone]] .& [x in P for x in df_plants[:,:GenID]]])
		#conv_units_in_zone = unique(df_plants.GenID[[x == zone_temp for x in df_plants[:,:Zone]] .& [x in P for x in df_plants[:,:GenID]]])
		if n in conv_nodes_in_zone
			conv_pmax_in_zone = sum(df_plants.Pmax[[x in conv_nodes_in_zone for x in df_plants[:,:OnBus]] .& [x in P for x in df_plants[:,:GenID]]])
			conv_pmax_at_node = sum(df_plants.Pmax[[x == n for x in df_plants[:,:OnBus]] .& [x in P for x in df_plants[:,:GenID]]])
			gsk_value_temp = conv_pmax_at_node/conv_pmax_in_zone
			gsk_temp[findfirst(N_FBMC .== n), findfirst(Z_FBMC .== zone_temp)] = gsk_value_temp
		end
	end
	return gsk_temp
end

gsk_pmax = get_gsk_pmax()
sum(gsk_pmax, dims=1)

function get_gsk_pmax_sub()
	P_sub = df_plants.GenID[[x in ["Hard Coal", "Gas/CCGT"] for x in df_plants[:,:Type]] .&
						    [x in Z_FBMC for x in df_plants[:,:Zone]]]
	gsk_temp = zeros(Float64, length(N_FBMC), length(Z_FBMC))
	for n in N_FBMC
		zone_temp = df_bus.Zone[df_bus[:,:BusID].==n][1]
		conv_nodes_in_zone = unique(df_plants.OnBus[[x == zone_temp for x in df_plants[:,:Zone]] .& [x in P_sub for x in df_plants[:,:GenID]]])
		#conv_units_in_zone = unique(df_plants.GenID[[x == zone_temp for x in df_plants[:,:Zone]] .& [x in P for x in df_plants[:,:GenID]]])
		if n in conv_nodes_in_zone
			conv_pmax_in_zone = sum(df_plants.Pmax[[x in conv_nodes_in_zone for x in df_plants[:,:OnBus]] .& [x in P_sub for x in df_plants[:,:GenID]]])
			conv_pmax_at_node = sum(df_plants.Pmax[[x == n for x in df_plants[:,:OnBus]] .& [x in P_sub for x in df_plants[:,:GenID]]])
			gsk_value_temp = conv_pmax_at_node/conv_pmax_in_zone
			gsk_temp[findfirst(N_FBMC .== n), findfirst(Z_FBMC .== zone_temp)] = gsk_value_temp
		end
	end
	return gsk_temp
end

gsk_pmax_sub = get_gsk_pmax_sub()
sum(gsk_pmax_sub, dims=1)

println("All GSKs built, with column sums")
println("1) GSK flat: ", round.(sum(gsk_flat,dims=1), digits=2))
println("2) GSK flat unit: ", round.(sum(gsk_flat_unit,dims=1), digits=2))
println("3) GSK pmax: ", round.(sum(gsk_pmax,dims=1), digits=2))
println("4) GSK pmax sub: ", round.(sum(gsk_pmax_sub,dims=1), digits=2))
