# Build zonal PTDFs and CNE selection
PTDF_Z = PTDF*gsk_cne

z2z_temp = zeros(Float64, length(L), Int(length(Z_FBMC)*(length(Z_FBMC)-1)/2))
for a in 1:1
	counter = 1
	for z in 1:(length(Z_FBMC)-1)
		for zz in (1+z):length(Z_FBMC)
			println("Zone 1: ", z, ", Zone 2: ", zz, ", Counter: ", counter)
			z2z_temp[:,counter] = PTDF_Z[:,z] - PTDF_Z[:,zz]
			counter = counter + 1
		end
	end
end

z2z_temp_abs = broadcast(abs, z2z_temp)
maximum_abs_z2z = reshape(findmax(z2z_temp_abs, dims=2)[1], :)

CNEC = L[findall(x->x>=cne_alpha, maximum_abs_z2z)]

if include_cb_lines==true
	cross_border_line=find_cross_border_lines()
	cross_border_line[findall(x-> !(x in CNEC), cross_border_line)]
	CNEC = unique(vcat(CNEC, cross_border_line))
end

# The order in CNEC does not necessarily follow order in L_in_CWE,
# therefore use CNEC = L[findall(...)]

CNEC = L[findall(x->x in CNEC, L)]

PTDF_CNEC = PTDF[findall(x->x in CNEC, L),:]
PTDF_Z_CNEC = PTDF_Z[findall(x->x in CNEC, L),:]

println("CNE selection: ",
		length(CNEC), " CNEs selected at alpha=", cne_alpha)
