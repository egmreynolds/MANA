#Functions used to summarise GWAS outputs
function calc_summaries(submodel)
    calc_tstats!(submodel)
    calc_pvalues_norm!(submodel)
end

function calc_tstats!(submodel)
    idx = findall(x -> x != -9, submodel.output_stats[:,1]) # all where MAF != -9
    if submodel.type == "additive"
        submodel.output_stats[idx,5] = submodel.output_stats[idx,3] ./ sqrt.(submodel.output_stats[idx,4])
    else
        submodel.output_stats[idx,7] = submodel.output_stats[idx,3] ./ sqrt.(submodel.output_stats[idx,4])
        submodel.output_stats[idx,8] = submodel.output_stats[idx,5] ./ sqrt.(submodel.output_stats[idx,6])
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
    end
end

function write_summaries(submodel)
    open(string(submodel.output_filename, ".part", submodel.id), "a") do io
        writedlm(io, [submodel.snpset.markerInfo submodel.output_stats], '\t'; header = false)
    end
end

function concatenate_summaries(model)
    fulloutput = Any["CHR" "BP" "SNP" "REF" "ALT" "MAF" "K" "BETA_A" "VAR_A" "BETA_D" "VAR_D" "T_A" "T_D" "P_A" "P_D"] # may need a check
    if model.type == "additive"
        fulloutput = Any["CHR" "BP" "SNP" "REF" "ALT" "MAF" "K" "BETA_A" "VAR_A" "T_A" "P_A"] # may
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
