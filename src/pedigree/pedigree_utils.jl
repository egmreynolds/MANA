# Functions to manage pedigree models, primarily utilise JWAS package - Utilities for pedigree_run.jl
function get_namesA(markerInfo)
     return markerInfo[:,3] .* "_" .* markerInfo[:,1] .* "_" .* markerInfo[:,2] .* "_A"
end

function make_genotype_table(genoprefix)
    geno = read_binary_genotypes(genoprefix)
    markerInfo = geno.markerInfo
    namesA = get_namesA(markerInfo)
    g_tbl_names = vcat(["anml_key"], namesA)
    g_tbl_names = map(x -> Symbol(x), g_tbl_names)
    return table(map(x -> parse(Int64,x), geno.IID), bedvals[geno.binarygenotypes][:]; names = g_tbl_names), markerInfo
end
#
function join_tables(y_tbl, g_tbl)
    return join(y_tbl, g_tbl; how = :inner, lkey = :anml_key, rkey = :anml_key)
end

function center_phenotype(tbl)
    meanY = mean(JuliaDB.select(tbl, :yd_trait))
    tbl = IndexedTables.transform(tbl, :yd_trait => :yd_trait => x -> x - meanY)
    return tbl
end
#
function generate_genotypeclass_columns(tbl, markerInfo)
    m = size(markerInfo, 1)
    namesA = get_namesA(markerInfo)
    namesG11 = markerInfo[:,3] .* "_" .* markerInfo[:,1] .* "_" .* markerInfo[:,2] .* "_G11"
    namesG12 = markerInfo[:,3] .* "_" .* markerInfo[:,1] .* "_" .* markerInfo[:,2] .* "_G12"
    namesG22 = markerInfo[:,3] .* "_" .* markerInfo[:,1] .* "_" .* markerInfo[:,2] .* "_G22"
    tbl = IndexedTables.transform(tbl, Symbol(namesG11[1]) => Symbol(namesA[1]) => x -> x == 0 ? 1 : 0)
    tbl = IndexedTables.transform(tbl, Symbol(namesG12[1]) => Symbol(namesA[1]) => x -> x == 1 ? 1 : 0)
    tbl = IndexedTables.transform(tbl, Symbol(namesG22[1]) => Symbol(namesA[1]) => x -> x == 2 ? 1 : 0)
    tbl = JuliaDB.select(tbl, (:anml_key, Symbol(namesG11[1]), Symbol(namesG12[1]), Symbol(namesG22[1])))
    return tbl
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
            chain = input[1][201:end]# remove burnin
            summary[i,3] = mean(chain)
            summary[i,4] = std(chain)
            summary[i,5] = calc_pvalue_norm(summary[i,3] ./ summary[i,4])
            summary[i,6] = gewekediag(chain)[2]
        end
        writedlm(basedir .* "/summary." .* pheno .* ".spliceQTL.txt", summary, '\t')
    end
end
