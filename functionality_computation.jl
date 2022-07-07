using ConScape
using Plots

thetas = (1.0, 0.1, 0.001)
scals = (224.3700107653066, 519.6098641681106, 5760.436073820083)
scals2 = (1, 1.5, 2, 3)

indir = "/output/"
outdir = "/output/"

for i in 1:10
    mov_prob, meta_p = ConScape.readasc(joinpath(indir, replace("ssf_sum_areaXXX.asc", "XXX" => i)))
    hab_qual, meta_q = ConScape.readasc(joinpath(indir, replace("rsf_sum_areaXXX.asc", "XXX" => i)))

    g = ConScape.Grid(size(mov_prob)...,
                        affinities=ConScape.graph_matrix_from_raster(mov_prob),
                        qualities=hab_qual,
                        costs=ConScape.mapnz(x -> -log(x), ConScape.graph_matrix_from_raster(mov_prob)));

    g_coarse = ConScape.Grid(size(mov_prob)...,
                        affinities=ConScape.graph_matrix_from_raster(mov_prob),
                        source_qualities=hab_qual,
                        target_qualities=ConScape.coarse_graining(g, 5),
                        costs=ConScape.mapnz(x -> -log(x), ConScape.graph_matrix_from_raster(mov_prob)));

    for j in 1:3
        h = ConScape.GridRSP(g_coarse, θ=thetas[j]);
		
		targetidx, targetnodes = ConScape._targetidx_and_nodes(h.g)
		qˢ = [h.g.source_qualities[i] for i in h.g.id_to_grid_coordinate_list]
		qᵗ = [h.g.target_qualities[i] for i in targetidx];
		EC = ConScape.expected_cost(h)

        for k in 1:4
            wparam = string(thetas[j], "_", scals2[k])
            wparam = replace(wparam, "." => "")
			
			func_vec = vec(sum(qˢ .* map(t -> iszero(t) ? t : exp(-t/scals[j]*scals2[k]), EC) .* qᵗ', dims=2)); #sum over cols

			func = fill(NaN, h.g.nrows, h.g.ncols)
			for (i,v) in enumerate(func_vec)
				func[h.g.id_to_grid_coordinate_list[i]] = v
			end
			
            file_name = replace("wrein_functionality_sum_areaXXX_paramsYYY.asc", "XXX" => i)
            file_name = replace(file_name, "YYY" => wparam)
            ConScape.writeasc(joinpath(outdir, file_name), func, meta_p)
        end
    end
end

thetas = ("lcp", "euclid")
scals = (170.59741434903637, 162.78820596099706)

for i in 1:10
    mov_prob, meta_p = ConScape.readasc(joinpath(indir, replace("ssf_sum_areaXXX.asc", "XXX" => i)))
    hab_qual, meta_q = ConScape.readasc(joinpath(indir, replace("rsf_sum_areaXXX.asc", "XXX" => i)))

    g = ConScape.Grid(size(mov_prob)...,
                        affinities=ConScape.graph_matrix_from_raster(mov_prob),
                        qualities=hab_qual,
                        costs=ConScape.mapnz(x -> -log(x), ConScape.graph_matrix_from_raster(mov_prob)));

    g_coarse = ConScape.Grid(size(mov_prob)...,
                        affinities=ConScape.graph_matrix_from_raster(mov_prob),
                        source_qualities=hab_qual,
                        target_qualities=ConScape.coarse_graining(g, 5),
                        costs=ConScape.mapnz(x -> -log(x), ConScape.graph_matrix_from_raster(mov_prob)));
						
	targetidx, targetnodes = ConScape._targetidx_and_nodes(g_coarse)
	qˢ = [g_coarse.source_qualities[i] for i in g_coarse.id_to_grid_coordinate_list]
	qᵗ = [g_coarse.target_qualities[i] for i in targetidx];


    for j in 1:2	
		if (j==1) EC = ConScape.least_cost_distance(g_coarse) end		
		if (j==2) EC = [hypot(xyᵢ[1] - xyⱼ[1], xyᵢ[2] - xyⱼ[2]) for xyᵢ in g_coarse.id_to_grid_coordinate_list, xyⱼ in targetidx] end

        for k in 1:4
            wparam = string(thetas[j], "_", scals2[k])
            wparam = replace(wparam, "." => "")
			
			func_vec = vec(sum(qˢ .* map(t -> iszero(t) ? t : exp(-t/scals[j]*scals2[k]), EC) .* qᵗ', dims=2)); #sum over cols

			func = fill(NaN, g_coarse.nrows, g_coarse.ncols)
			for (i,v) in enumerate(func_vec)
				func[g_coarse.id_to_grid_coordinate_list[i]] = v
			end
			
            file_name = replace("wrein_functionality_sum_areaXXX_paramsYYY.asc", "XXX" => i)
            file_name = replace(file_name, "YYY" => wparam)
            ConScape.writeasc(joinpath(outdir, file_name), func, meta_p)
        end
    end
end

