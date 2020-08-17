using MANA
using JuliaDB, JWAS, DataFrames

function main()
	dir="~/MANA/example/"
        snpfilename=dir .* "plink.example"
        phenofilename="phen.example2.txt"
        pedfilename="ped.example2.csv"
        repeated_measures=false
        outputfilename="out.example2.txt"
        Rvar=0.7 #total, not a proportion
        Gvar=0.3 #total, not a proportion
	Pvar=0
        y_tbl = loadtable(phenofilename; spacedelim = true)
        pedigree = get_pedigree(pedfilename)
        g_tbl, markerInfo = make_genotype_table(snpfilename)
        g_tbl = generate_genotypeclass_columns(g_tbl, markerInfo)
        tbl = join(y_tbl, g_tbl; lkey=:anml_key, rkey=:anml_key)
        tbl = MANA.center_phenotype(tbl)
        covariates = []
        if repeated_measures
                tbl = IndexedTables.transform(tbl, :anml_keyx => :anml_key)
                df = DataFrame(tbl)
                run_pedigree_model_multirecord(df, markerInfo, pedigree, Rvar, Gvar, Pvar; chainlength = 50000, outputfilename = outputfilename, covariates = covariates)
        else
                df = DataFrame(tbl)
                run_pedigree_model(df, markerInfo, pedigree, Rvar, Gvar; chainlength = 50000, outputfilename = outputfilename, covariates = covariates)
        end
end

main()

