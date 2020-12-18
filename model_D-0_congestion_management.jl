cd("./results/D-1_market_coupling")

d_1_curt = convert(Array{Float64,2}, CSV.read(string("df_d_1_curt.csv")))
d_1_gen = convert(Array{Float64,2}, CSV.read(string("df_d_1_gen.csv")))
d_1_gen_costs = convert(Array{Float64,2}, CSV.read(string("df_d_1_gen_costs.csv")))
d_1_nodal_price= convert(Array{Float64,2}, CSV.read(string("df_d_1_nodal_price.csv")))

cd("..")
cd("..")

# Create empty matrices to store values
d_0_curt = zeros(Float64, length(T), length(N))
d_0_rd_pos = zeros(Float64, length(T), length(P_RD))
d_0_rd_neg = zeros(Float64, length(T), length(P_RD))
d_0_nod_inj = zeros(Float64, length(T), length(N))
d_0_line_f = zeros(Float64, length(T), length(L))
d_0_delta = zeros(Float64, length(T), length(N))
d_0_gen_costs_rd_pos = zeros(Float64, length(T), length(Z))
d_0_gen_costs_rd_neg = zeros(Float64, length(T), length(Z))
d_0_redispatch_costs = zeros(Float64, length(T), length(Z))
d_0_curt_costs = zeros(Float64, length(T), length(Z))

days_foresight = 1
hours_per_horizon = 4*168

max_mc = find_maximum_mc()
cost_curt = ceil(max_mc)

for horizon in 1:ceil(Int, length(T)/hours_per_horizon)

    println("Horizon: ", horizon, "/", ceil(Int, length(T)/hours_per_horizon))
	Tsub = ((horizon-1)*hours_per_horizon+1):min((horizon*hours_per_horizon+(days_foresight-1)*24), length(T))

	m = Model(Gurobi.Optimizer)

	@variable(m, -d_1_curt[t,findfirst(N .== n)] <= CURT_RD[t in Tsub, n in N] <= get_renew(t,n) - d_1_curt[t,findfirst(N .== n)])
	@variable(m, DELTA[t in Tsub, n in N])
	@variable(m, NOD_INJ[t in Tsub, n in N])
	@variable(m, LINE_F[t in Tsub, l in L])
	@variable(m, 0 <= RD_POS[t in Tsub, p in P_RD] <= get_gen_up(p) - d_1_gen[t,findfirst(P .== p)])
	@variable(m, 0 <= RD_NEG[t in Tsub, p in P_RD] <= d_1_gen[t,findfirst(P .== p)])
	@variable(m, GEN_COSTS_RD_POS[t in Tsub, z in Z])
	@variable(m, GEN_COSTS_RD_NEG[t in Tsub, z in Z])
	@variable(m, REDISPATCH_COSTS[t in Tsub, z in Z])
	@variable(m, CURT_COSTS_RD[t in Tsub, z in Z])
	@objective(m, Min,
	sum(RD_POS[t,p]*(100+get_mc(p)) for t in Tsub for p in P_RD) +
	sum(RD_NEG[t,p]*(100+cost_curt-get_mc(p)) for t in Tsub for p in P_RD) +
	sum(CURT_RD[t,n]*(100+cost_curt) for t in Tsub for n in N))

	println("Variables done.")

    @constraint(m, costs_gen_rd_pos[t=Tsub, z=Z],
    sum(sum(RD_POS[t,p]*get_mc(p) for p in p_rd_at_n[n]) for n in n_in_z[z])
    == GEN_COSTS_RD_POS[t,z])
    println("Built constraints costs_gen_rd_pos.")

	@constraint(m, costs_gen_rd_neg[t=Tsub, z=Z],
    sum(sum(RD_NEG[t,p]*get_mc(p) for p in p_rd_at_n[n]) for n in n_in_z[z])
    == GEN_COSTS_RD_NEG[t,z])
    println("Built constraints costs_gen_rd_neg.")

    @constraint(m, costs_redispatch[t=Tsub, z=Z],
    REDISPATCH_COSTS[t,z] == GEN_COSTS_RD_POS[t,z] - GEN_COSTS_RD_NEG[t,z])
    println("Built constraints costs_redispatch.")

    @constraint(m, costs_curt_rd[t=Tsub, z=Z],
    sum(CURT_RD[t,n]*cost_curt for n in n_in_z[z])
    == CURT_COSTS_RD[t,z])
    println("Built constraints costs_curt_rd.")

	#n=59
	#t=1
	#p=p_at_n[n]
	#d_1_gen[1:24,findfirst(P .== p)]
	#get_line_cap(1)

    @constraint(m, nodal_balance[t=Tsub, n=N],
    sum(d_1_gen[t,findfirst(P .== p)] for p in p_at_n[n])
    + sum(RD_POS[t,p] for p in p_rd_at_n[n])
    - sum(RD_NEG[t,p] for p in p_rd_at_n[n])
    - NOD_INJ[t,n]
    + get_renew(t,n)
    - d_1_curt[t,findfirst(N .== n)]
    - CURT_RD[t,n]
    ==
    get_dem(t,n))
    println("Built constraints nodal_balance.")

	@constraint(m, nodal_injection[t=Tsub, n=N],
	NOD_INJ[t,n] == sum(B_mat[n,nn]*DELTA[t,nn] for nn in N))
	println("Built constraints nodal_injection.")

	@constraint(m, line_flow[t=Tsub, l=L],
	LINE_F[t,l] == sum(H_mat[l,n]*DELTA[t,n] for n in N))
	println("Built constraints line_flow.")

    @constraint(m, line_cap_pos[t=Tsub, l=L], LINE_F[t,l] <= get_line_cap(l))
    println("Built constraints line_cap_pos.")

	@constraint(m, line_cap_neg[t=Tsub, l=L], -LINE_F[t,l] <= get_line_cap(l))
    println("Built constraints line_cap_neg.")

	for t in Tsub
		JuMP.fix(DELTA[t,68], 0)
	end
	println("Built constraints FIX SLACK NODE.")

    println("Constraints done.")

    status = optimize!(m)

    d_0_curt[Tsub,:] = JuMP.value.(CURT_RD[:,:])
    d_0_rd_pos[Tsub,:] = JuMP.value.(RD_POS[:,:])
    d_0_rd_neg[Tsub,:] = JuMP.value.(RD_NEG[:,:])
    d_0_gen_costs_rd_pos[Tsub,:] = JuMP.value.(GEN_COSTS_RD_POS[:,:])
    d_0_gen_costs_rd_neg[Tsub,:] = JuMP.value.(GEN_COSTS_RD_NEG[:,:])
    d_0_redispatch_costs[Tsub,:] = JuMP.value.(REDISPATCH_COSTS[:,:])
    d_0_curt_costs[Tsub,:] = JuMP.value.(CURT_COSTS_RD[:,:])
    d_0_nod_inj[Tsub,:] = JuMP.value.(NOD_INJ[:,:])
    d_0_line_f[Tsub,:] = JuMP.value.(LINE_F[:,:])
    d_0_delta[Tsub,:] = JuMP.value.(DELTA[:,:])

    m = nothing
end

include("export_D-0_results.jl")
