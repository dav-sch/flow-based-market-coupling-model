d_2_curt = zeros(Float64, length(T), length(N))
d_2_delta = zeros(Float64, length(T), length(N))
d_2_nod_inj = zeros(Float64, length(T), length(N))
d_2_line_f = zeros(Float64, length(T), length(L))
d_2_gen = zeros(Float64, length(T), length(P))
d_2_gen_costs = zeros(Float64, length(T), length(Z))
d_2_curt_costs = zeros(Float64, length(T), length(Z))
d_2_nodal_price = zeros(Float64, length(T), length(N))
d_2_np = zeros(Float64, length(T), length(Z_FBMC))
d_2_export = zeros(Float64, length(T), length(Z_not_in_FBMC))

max_mc = find_maximum_mc()
cost_curt = ceil(max_mc)
hours_per_horizon = 4*168
days_foresight = 1

for horizon in 1:ceil(Int, length(T)/hours_per_horizon)

	println("Horizon: ", horizon, "/", ceil(Int, length(T)/hours_per_horizon))
	Tsub = ((horizon-1)*hours_per_horizon+1):min((horizon*hours_per_horizon+(days_foresight-1)*24), length(T))

	m = Model(Gurobi.Optimizer)

	@variable(m, 0 <= CURT[t in Tsub, n in N] <= get_renew(t,n))
	@variable(m, DELTA[t in Tsub, n in N])
	@variable(m, NOD_INJ[t in Tsub, n in N])
	@variable(m, LINE_F[t in Tsub, l in L])
	@variable(m, 0 <= GEN[t in Tsub, p in P] <= get_gen_up(p))
	@variable(m, GEN_COSTS[t in Tsub, z in Z])
	@variable(m, CURT_COSTS[t in Tsub, z in Z])
	@variable(m, NP[t in Tsub, z in Z_FBMC])
	@variable(m, EXPORT[t in Tsub, z in Z_not_in_FBMC])
	@objective(m, Min,
		sum(GEN[t,p]*get_mc(p) for t in Tsub for p in P) +
		sum(CURT[t,n]*cost_curt for t in Tsub for n in N))
	println("Variables done.")

	@constraint(m, costs_gen[t=Tsub, z=Z], sum(GEN[t,p]*get_mc(p) for p in p_in_z[z]) == GEN_COSTS[t,z])
	println("Built constraints costs_gen.")

	@constraint(m, costs_curt[t=Tsub, z=Z], sum(CURT[t,n]*cost_curt for n in n_in_z[z]) == CURT_COSTS[t,z])
	println("Built constraints costs_curt.")

	@constraint(m, nodal_balance[t=Tsub, n=N],
	sum(GEN[t,p] for p in p_at_n[n]) + get_renew(t,n)
	- NOD_INJ[t,n] - CURT[t,n]
	==
	get_dem(t,n))
	println("Built constraints nodal_balance.")

	@constraint(m, export_balance_abroad[t=Tsub, z=Z_not_in_FBMC],
	EXPORT[t,z] == sum(NOD_INJ[t,n] for n in n_in_z[z]))
	println("Export balance outside of FBMC.")

	@constraint(m, net_positions[t=Tsub, z=Z_FBMC],
	NP[t,z] == sum(NOD_INJ[t,n] for n in n_in_z[z]) +
		       sum(EXPORT[t,zz] for zz in ifelse(z in Z_FBMC, z_to_z[z][findall(z->!(z in Z_FBMC), z_to_z[z])], z_to_z[z])))
	println("Net positions within FBMC.")

	if np_within_fbmc_0 == true
		@constraint(m, no_trade_base_case_FBMC[t=Tsub, z=Z_FBMC], NP[t,z] == 0)
		println("No trade condition (within FBMC).")
	end

	if exports_abroad_0 == true
		@constraint(m, no_trade_base_case_abroad[t=Tsub, z=Z_not_in_FBMC], EXPORT[t,z] == 0)
		println("No trade condition (within FBMC).")
	end

	@constraint(m, trade_limit_base_case_pos[t=Tsub, z=Z_not_in_FBMC], EXPORT[t,z] <= max_ntc)
	println("Trade limit pos. (outside FBMC).")

	@constraint(m, trade_limit_base_case_neg[t=Tsub, z=Z_not_in_FBMC], EXPORT[t,z] >= -max_ntc)
	println("Trade limit neg. (outside FBMC).")

	@constraint(m, nodal_injection[t=Tsub, n=N],
	NOD_INJ[t,n] == sum(B_mat[n,nn]*DELTA[t,nn] for nn in N))
	println("Built constraints nodal_injection.")

	@constraint(m, line_flow[t=Tsub, l=L],
	LINE_F[t,l] == sum(H_mat[l,n]*DELTA[t,n] for n in N))
	println("Built constraints line_flow.")

	@constraint(m, line_cap_pos[t=Tsub, l=L], LINE_F[t,l] <= get_line_cap(l))
	println("Built constraints line_cap_pos.")

	@constraint(m, line_cap_neg[t=Tsub, l=L], LINE_F[t,l] >= -get_line_cap(l))
	println("Built constraints line_cap_neg.")

	@constraint(m, line_cap_pos_cnec[t=Tsub, l=CNEC],
	LINE_F[t,l] <= get_line_cap(l)*(1-frm))
	println("Built constraints line_cap_pos for CNECs.")

	@constraint(m, line_cap_neg_cnec[t=Tsub, l=CNEC],
	LINE_F[t,l] >= -get_line_cap(l)*(1-frm))
	println("Built constraints line_cap_neg for CNECs.")

	for t in Tsub
		JuMP.fix(DELTA[t,68], 0)
	end
	println("Built constraints FIX SLACK NODE.")

	println("Constraints done.")

	status = optimize!(m)

	d_2_curt[Tsub,:] = JuMP.value.(CURT[:,:])
	d_2_delta[Tsub,:] = JuMP.value.(DELTA[:,:])
	d_2_nod_inj[Tsub,:] = JuMP.value.(NOD_INJ[:,:])
	d_2_line_f[Tsub,:] = JuMP.value.(LINE_F[:,:])
	d_2_gen[Tsub,:] = JuMP.value.(GEN[:,:])
	d_2_gen_costs[Tsub,:] = JuMP.value.(GEN_COSTS[:,:])
	d_2_curt_costs[Tsub,:] = JuMP.value.(CURT_COSTS[:,:])
	d_2_np[Tsub,:] = JuMP.value.(NP[:,:])
	d_2_export[Tsub,:] = JuMP.value.(EXPORT[:,:])
	d_2_nodal_price[Tsub,:] = dual.(nodal_balance[:,:])

	m = nothing
end

include("export_D-2_results.jl")
