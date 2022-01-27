# Functions to manage pedigree models, primarily utilise JWAS package - Utilities for pedigree_run.jl
function get_namesA(markerInfo)
     return markerInfo[:,3] .* "_" .* markerInfo[:,1] .* "_" .* markerInfo[:,2] .* "_A"
end

function make_genotype_table(genoprefix)
    geno = read_binary_genotypes(genoprefix)
    df = DataFrame(anml_key = map(x -> parse(Int64,x), geno.IID), G = bedvals[geno.binarygenotypes[:,1]])
    return df, geno.markerInfo
end
#
function join_tables(y_tbl, g_tbl)
    return join(y_tbl, g_tbl; how = :inner, lkey = :anml_key, rkey = :anml_key)
end

function center_phenotype!(df)
     df.yd_trait = df.yd_trait .- mean(df.yd_trait)
end
#
function generate_genotypeclass_columns(df, markerInfo)
     df.G11 = map(x -> x == 0 ? 1 : 0, df.G)
     df.G12 = map(x -> x == 1 ? 1 : 0, df.G)
     df.G22 = map(x -> x == 2 ? 1 : 0, df.G)
     select!(df, [:anml_key, :G11, :G12, :G22])
     return df
end

function calc_pvalue_norm(t)
    dist = Normal()
    return 2*ccdf.(dist, abs.(t))
end

function summarise(pheno_header, basedir, markerlist; burnin=200)
    m = size(markerlist, 1)
    for pheno in pheno_header
        dir = basedir .* "/" .* pheno
        summary = [repeat([pheno],m) markerlist zeros(Float64, m, 4)]
        for i in 1:m
            file = dir .* "/" .* markerlist[i] .* "_1:" .* markerlist[i] .* ".txt"
            input = readdlm(file, '\t'; header = true)
            chain = input[1][(burnin+1):end] # remove burnin
            summary[i,3] = mean(chain)
            summary[i,4] = std(chain)
            summary[i,5] = calc_pvalue_norm(summary[i,3] ./ summary[i,4])
            summary[i,6] = gewekediag(chain)[2]
        end
        writedlm(basedir .* "/summary." .* pheno .* ".spliceQTL.txt", summary, '\t')
    end
end
