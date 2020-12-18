cd("./results/D-2_base_case")

d_2_line_f = convert(Array{Float64,2}, CSV.read(string("df_d_2_line_f.csv")))
d_2_np = convert(Array{Float64,2}, CSV.read(string("df_d_2_np.csv")))
d_2_gen = convert(Array{Float64,2}, CSV.read(string("df_d_2_gen.csv")))

cd("..")
cd("..")

PTDF_Z = PTDF*gsk_mc
PTDF_Z_CNEC = PTDF_Z[findall(x->x in CNEC, L),:]

# Create empty matrices to store values
d_1_curt = zeros(Float64, length(T), length(N))
d_1_gen = zeros(Float64, length(T), length(P))
d_1_np = zeros(Float64, length(T), length(Z_FBMC))
d_1_dump_dem = zeros(Float64, length(T), length(Z))
d_1_gen_costs = zeros(Float64, length(T), length(Z))
d_1_curt_costs = zeros(Float64, length(T), length(Z))
d_1_nodal_price = zeros(Float64, length(T), length(Z_FBMC))
d_1_nodal_price_abroad = zeros(Float64, length(T), length(Z_not_in_FBMC))

hours_per_horizon = 4*168
days_foresight = 1

for horizon in 1:ceil(Int, length(T)/hours_per_horizon)

	println("Horizon: ", horizon, "/", ceil(Int, length(T)/hours_per_horizon))
	Tsub = ((horizon-1)*hours_per_horizon+1):min((horizon*hours_per_horizon+(days_foresight-1)*24), length(T))

	m = Model(Gurobi.Optimizer)

	@variable(m, 0 <= CURT[t in Tsub, n in N] <= get_renew(t,n))
	@variable(m, NP[t in Tsub, z in Z_FBMC])
	@variable(m, 0 <= GEN[t in Tsub, p in P] <= get_gen_up(p))
	@variable(m, 0 <= EXPORT[t in Tsub, z in Z, zz in z_to_z[z]] <= max_ntc)
	@variable(m, GEN_COSTS[t in Tsub, z in Z])
	@variable(m, CURT_COSTS[t in Tsub, z in Z])
	@objective(m, Min,
		sum(GEN[t,p]*get_mc(p) for t in Tsub for p in P)
		+ sum(CURT[t,n]*cost_curt_mc for t in Tsub for n in N)
		)
	println("Variables done.")

	@constraint(m, costs_gen[t=Tsub, z=Z],
	sum(GEN[t,p]*get_mc(p) for p in p_in_z[z]) == GEN_COSTS[t,z])
	println("Built constraints costs_gen.")

	@constraint(m, costs_curt[t=Tsub, z=Z],
	sum(CURT[t,n]*cost_curt_mc for n in n_in_z[z]) == CURT_COSTS[t,z])
	println("Built constraints costs_curt.")

	@constraint(m, zonal_balance[t=Tsub, z=Z_FBMC],
		sum(GEN[t,p] for p in p_in_z[z])
		+ sum(get_renew(t,n) for n in n_in_z[z])
		- sum(CURT[t,n] for n in n_in_z[z])
		- sum(EXPORT[t,z,zz] for zz in ifelse(z in Z_FBMC, z_to_z[z][findall(z->!(z in Z_FBMC), z_to_z[z])], z_to_z[z]))
		+ sum(EXPORT[t,zz,z] for zz in ifelse(z in Z_FBMC, z_to_z[z][findall(z->!(z in Z_FBMC), z_to_z[z])], z_to_z[z]))
		- NP[t,z]
		==
		sum(get_dem(t,n) for n in n_in_z[z])
		)
	println("Built constraints zonal_balance.")

	@constraint(m, zonal_balance_abroad[t=Tsub, z=Z_not_in_FBMC],
		sum(GEN[t,p] for p in p_in_z[z])
		+ sum(get_renew(t,n) for n in n_in_z[z])
		- sum(CURT[t,n] for n in n_in_z[z])
		- sum(EXPORT[t,z,zz] for zz in ifelse(z in Z_FBMC, z_to_z[z][findall(z->!(z in Z_FBMC), z_to_z[z])], z_to_z[z]))
		+ sum(EXPORT[t,zz,z] for zz in ifelse(z in Z_FBMC, z_to_z[z][findall(z->!(z in Z_FBMC), z_to_z[z])], z_to_z[z]))
		==
		sum(get_dem(t,n) for n in n_in_z[z])
		#- DUMP_DEM[t,z]
		)
	println("Built constraints zonal_balance abroad.")

	@constraint(m, net_position_fbmc[t=Tsub],
	sum(NP[t,z] for z in Z_FBMC) == 0)
	println("Net positions sum-zero inside of FBMC.")

	@constraint(
	m, flow_on_cnes_pos[t=Tsub, j=CNEC],
	sum(PTDF_Z_CNEC[findfirst(CNEC .== j),findfirst(Z_FBMC .== z_fb)]*
		(NP[t,z_fb]-d_2_np[t,findfirst(Z .== z_fb)]) for z_fb in Z_FBMC)
	<= get_line_cap(j)*(1-frm)-d_2_line_f[t,findfirst(L .== j)]
	)
	println("Flows on CNECs (pos).")

	@constraint(
	m, flow_on_cnes_neg[t=Tsub, j=CNEC],
	sum(PTDF_Z_CNEC[findfirst(CNEC .== j),findfirst(Z_FBMC .== z_fb)]*
		(NP[t,z_fb]-d_2_np[t,findfirst(Z .== z_fb)]) for z_fb in Z_FBMC)
	>= -get_line_cap(j)*(1-frm)-d_2_line_f[t,findfirst(L .== j)]
	)
	println("Flows on CNECs (neg).")

	println("Constraints done.")

	status = optimize!(m)

	d_1_curt[Tsub,:] = JuMP.value.(CURT[:,:])
	d_1_np[Tsub,:] = JuMP.value.(NP[:,:])
	d_1_gen[Tsub,:] = JuMP.value.(GEN[:,:])
	d_1_gen_costs[Tsub,:] = JuMP.value.(GEN_COSTS[:,:])
	d_1_curt_costs[Tsub,:] = JuMP.value.(CURT_COSTS[:,:])
	d_1_nodal_price[Tsub,:] = dual.(zonal_balance[:,:])
	d_1_nodal_price_abroad[Tsub,:] = dual.(zonal_balance_abroad[:,:])

	m = nothing
end

include("export_D-1_results.jl")
