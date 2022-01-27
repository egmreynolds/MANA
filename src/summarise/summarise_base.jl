#Functions used to summarise GWAS outputs
function calc_summaries(submodel)
    calc_zstats!(submodel)
    calc_pvalues_norm!(submodel)
end

function calc_zstats!(submodel)
    idx = findall(x -> x != -9, submodel.output_stats[:,1]) # all where MAF != -9
    if submodel.type == "additive"
        submodel.output_stats[idx,5] = submodel.output_stats[idx,3] ./ sqrt.(submodel.output_stats[idx,4])
    else
        submodel.output_stats[idx,7] = submodel.output_stats[idx,3] ./ sqrt.(submodel.output_stats[idx,4])
        submodel.output_stats[idx,8] = submodel.output_stats[idx,5] ./ sqrt.(submodel.output_stats[idx,6])
	submodel.output_stats[idx,14] = submodel.output_stats[idx,12] ./ sqrt.(submodel.output_stats[idx,13])
        submodel.output_stats[idx,20] = submodel.output_stats[idx,16] ./ sqrt.(submodel.output_stats[idx,17])
        submodel.output_stats[idx,21] = submodel.output_stats[idx,18] ./ sqrt.(submodel.output_stats[idx,19])
    end
end

function calc_pvalues_norm!(submodel)
    dist = Normal()
    idx = findall(x -> x != -9, submodel.output_stats[:,1])
    if submodel.type == "additive"
        submodel.output_stats[idx,6] = 2*ccdf.(dist, abs.(submodel.output_stats[idx,5]))
    else
        submodel.output_stats[idx,9] = 2*ccdf.(dist, abs.(submodel.output_stats[idx,7]))
        submodel.output_stats[idx,10] = 2*ccdf.(dist, abs.(submodel.output_stats[idx,8]))
	submodel.output_stats[idx,15] = 2*ccdf.(dist, abs.(submodel.output_stats[idx,14]))
        submodel.output_stats[idx,22] = 2*ccdf.(dist, abs.(submodel.output_stats[idx,20]))
        submodel.output_stats[idx,23] = 2*ccdf.(dist, abs.(submodel.output_stats[idx,21]))
    end
end

function write_summaries(submodel)
    open(string(submodel.output_filename, ".part", submodel.id), "a") do io
        writedlm(io, [submodel.snpset.markerInfo submodel.output_stats], '\t'; header = false)
    end
end

function concatenate_summaries(model)
    fulloutput = Any["CHR" "BP" "SNP" "REF" "ALT" "MAF" "K" "BETA_A" "VAR_A" "BETA_D" "VAR_D" "Z_A" "Z_D" "P_A" "P_D" "BETA_I" "BETA_ALPHA" "VAR_ALPHA" "Z_ALPHA" "P_ALPHA" "BETA_HET" "VAR_HET" "BETA_HOM" "VAR_HOM" "Z_HET" "Z_HOM" "P_HET" "P_HOM"] 
    if model.type == "additive"
        fulloutput = Any["CHR" "BP" "SNP" "REF" "ALT" "MAF" "K" "BETA_A" "VAR_A" "Z_A" "P_A"] # may
    end
    for future in model.futures
        wait(future)
    end
    for i in 1:model.submodelcount
        filename = string(model.outputfilename, ".part", (i-1))
	if filesize(filename) > 0
	        fulloutput = vcat(fulloutput, readdlm(filename, '\t', Any, '\n'; header = false))
	end
        run(`rm $filename`)
    end
    open(string(model.outputfilename, ".txt"), "a") do io
        writedlm(io, fulloutput, '\t'; header = false)
    end
end
