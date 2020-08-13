# Functions to run pedigree models
function run_pedigree_model_multirecord(tbl, markerInfo, pedigree, Rvar, Gvar, Pvar; chainlength = 10000, outputfilename = "MCMC_samples", covariates = [])
    namesG11 = markerInfo[:,3] .* "_" .* markerInfo[:,1] .* "_" .* markerInfo[:,2] .* "_G11"
    namesG12 = markerInfo[:,3] .* "_" .* markerInfo[:,1] .* "_" .* markerInfo[:,2] .* "_G12"
    namesG22 = markerInfo[:,3] .* "_" .* markerInfo[:,1] .* "_" .* markerInfo[:,2] .* "_G22"
    covars = [namesG11[1]; namesG12[1]; namesG22[1]; covariates]
    printcovars = [covars; "anml_keyx"]
    eqcovars = [covars; "anml_key"; "anml_keyx"]
    model_equation = "yd_trait ="
    model_terms = join(eqcovars, " + ")
    model_equation = join([model_equation model_terms], " ")
    print("MODEL EQUATION:  ")
    println(model_equation)
    model = build_model(model_equation, Rvar)
    for term in covars
        set_covariate(model, term)
    end
    set_random(model,"anml_key",pedigree,Gvar)
    set_random(model, "anml_keyx", Pvar)
    for term in printcovars
        outputMCMCsamples(model, term)
    end
    out = runMCMC(model, tbl, output_samples_frequency = 10, chain_length = chainlength, output_samples_file = outputfilename)
end

function run_pedigree_model(tbl, markerInfo, pedigree, Rvar, Gvar; chainlength = 10000, outputfilename = "MCMC_samples", covariates = [])
    namesG11 = markerInfo[:,3] .* "_" .* markerInfo[:,1] .* "_" .* markerInfo[:,2] .* "_G11"
    namesG12 = markerInfo[:,3] .* "_" .* markerInfo[:,1] .* "_" .* markerInfo[:,2] .* "_G12"
    namesG22 = markerInfo[:,3] .* "_" .* markerInfo[:,1] .* "_" .* markerInfo[:,2] .* "_G22"
    covars = [namesG11[1]; namesG12[1]; namesG22[1]; covariates]
    printcovars = covars
    eqcovars = [covars; "anml_key"]
    model_equation = "yd_trait ="
    model_terms = join(eqcovars, " + ")
    model_equation = join([model_equation model_terms], " ")
    print("MODEL EQUATION:  ")
    println(model_equation)
    model = build_model(model_equation, Rvar)
    for term in covars
        set_covariate(model, term)
    end
    set_random(model,"anml_key",pedigree,Gvar)
    for term in printcovars
        outputMCMCsamples(model, term)
    end
    out = runMCMC(model, tbl, output_samples_frequency = 10, chain_length = chainlength, output_samples_file = outputfilename)
end

function run_pedigree_model_eQTL(tbl, pedigree, Rvar, Gvar; chainlength = 10000, freq = 10, outputfilename = "eQTL.MCMC_samples")
    pheno = String(names(tbl)[2])
    snp = String(names(tbl)[3])
    model_equation = pheno .* " = intercept + " .* snp .* " + anml_key"
    model = build_model(model_equation, Rvar)
    set_covariate(model, snp)
    set_random(model,"anml_key",pedigree, Gvar)
    outputMCMCsamples(model, snp)
    out = runMCMC(model, tbl, output_samples_frequency = freq, chain_length = chainlength, output_samples_file = outputfilename)
end

function run_spliceQTL(snp_data, pheno_data, pedigree, Rvar, Gvar; chainlength = 10000, freq = 10, outdir = "MCMC_samples")
    k = size(pheno_data.header, 1)
    m = size(snp_data.header, 1)
    for i in 1:k
        pheno_symbol = Symbol(pheno_data.header[i])
        phenoi = table(pheno_data.iid, pheno_data.data[:,i], names = [:anml_key, pheno_symbol], pkey=:anml_key)
        outputdir = outdir .* "/" .* pheno_data.header[i]
        run(`mkdir -p $outputdir`)
        for j in 1:m
            outputfilename = outputdir .* "/" .* snp_data.header[j]
            snp_symbol = Symbol(snp_data.header[j])
            if any(snp_data.data[:,j] .== "NA") # change NAs to mean(genotype)
                na_index = findall(snp_data.data[:,j] .== "NA")
                snp_data.data[na_index,j] = repeat([missing], size(na_index,1))
                snp_data.data[na_index,j] = repeat([mean(skipmissing(snp_data.data[:,j]))], size(na_index,1))
            end
            snpj = table(snp_data.iid, snp_data.data[:,j], names = [:anml_key, snp_symbol], pkey = :anml_key)
            tbl = join_tables(phenoi, snpj)
            meanA = mean(IndexedTables.select(tbl, snp_symbol))
            tbl = IndexedTables.transform(tbl, snp_symbol => snp_symbol => x -> x - meanA)
            df = DataFrame(tbl)
            run_pedigree_model_eQTL(df, pedigree, Rvar, Gvar; chainlength = chainlength, freq = freq, outputfilename = outputfilename)
        end
    end
end