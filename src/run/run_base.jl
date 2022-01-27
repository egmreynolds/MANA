# Functions for running single-locus non-additive and additive gwas models
function run_model(model)
    m1 = 1
    m2 = model.setsize
    while model.notcomplete
        if m2 >= model.m
            m2 = model.m
            model.notcomplete = false
        end
        curridx = 1
        snpset = get_bed_values(model.snps, m1, m2)
        snpset_maxbp = parse(Float64, snpset.markerInfo[end,2])
        makesubmodels = true
        while makesubmodels
            if snpset_maxbp <= model.layerlimit
                snpset = get_bed_values(model.snps, m1 + curridx - 1, m2)
                job = @spawn run_association(SubModel(size(snpset.markerInfo, 1), model.n, model.k, snpset, model.phenotypes, model.phenotypesfilename, nothing, nothing, nothing, nothing, model.outputfilename, model.printchain, model.submodelcount, model.weights, model.type))
                push!(model.futures, job)
                sleep(1)
                model.submodelcount += 1
                makesubmodels = false
            else
                bpidx_less = findall(x -> parse(Float64, x) <= model.layerlimit, snpset.markerInfo[curridx:end,2])
                try
                    snpset_less = get_bed_values(model.snps, m1 + curridx - 1, m1 + bpidx_less[end] - 1)
                    job = @spawn run_association(SubModel(size(snpset_less.markerInfo,1), model.n, model.k, snpset_less, model.phenotypes, model.phenotypesfilename, nothing, nothing, nothing, nothing, model.outputfilename, model.printchain, model.submodelcount, model.weights, model.type))
                    push!(model.futures, job)
                    sleep(1)
                    curridx = 1 + bpidx_less[end]
                    model.submodelcount += 1
                catch error
                    if !isa(error, BoundsError)
                        println("ERROR: Unknown error")
                        println(error)
                    else
                        println("BoundsError: Switching to next phenotypes interval")
                        println(model.submodelcount)
                        next_phenotype_interval!(model)
                        continue
                    end
                end

                if curridx <= (m2-m1+1)
                    # switch phenotype
                    next_phenotype_interval!(model)
                else
                    # you reach the end of setsize and still in same pheno, so move on to next set.
                    makesubmodels = false
                end
            end
        end
        m1 += model.setsize
        m2 += model.setsize
    end
end

function run_association(submodel)
    inititate_submodel!(submodel)
    if submodel.type == "additive"
        for i in 1:submodel.m
            testsnp = select_snp(submodel, i)
            calc_additive_only_stats!(submodel, testsnp, i)
        end
    else
        for i in 1:submodel.m
            testsnp = select_snp(submodel, i)
            calc_dominance_stats!(submodel, testsnp, i)
        end

    end
    calc_summaries(submodel)
    write_summaries(submodel)
end

function calc_dominance_stats!(submodel, testsnp, i)
    cellsums, cellcounts, zero_cells = calc_betas(testsnp.genotypes[:], submodel.phenotype.trait, submodel.n, submodel.weights.rinverse)
    if zero_cells
        submodel.output_stats[i,:] = fill(-9, 23)
        return
    end
    invcellcounts = 1 ./ cellcounts
    betahat = cellsums .* invcellcounts
    mse = diag(submodel.yty .- betahat'cellsums) ./ (submodel.n-3)
    sigma2 = invcellcounts*mse'
    sigma = sqrt.(sigma2)
    z = randn(submodel.rng, Float64, (3,submodel.k))
    beta_plaus = (z .* sigma) .+ betahat
    if submodel.printchain
	writedlm(submodel.output_filename .* "." .* submodel.snpset.markerInfo[i,3] .* ".tsv", beta_plaus', '\t')
    end
    
    maf = sum(Array{Int64,1}(testsnp.genotypes[:]))/(2*submodel.n)
    K = Array{Float64,2}([1.0 0.0 0.0
    -0.5 0.0 0.5
    -0.5 1.0 -0.5
    maf-1 1-2*maf maf
    -1.0 1.0 0.0
    -1.0 0.0 1.0])

    beta = K*beta_plaus		

    aview = view(beta, 2, :)
    dview = view(beta, 3, :)
    #
    iview = view(beta, 1, :)
    alphaview = view(beta, 4, :)
    hetview = view(beta, 5, :)
    homview = view(beta, 6, :)
    mean_i = mean(iview)
    mean_alpha = mean(alphaview)
    var_alpha = var(alphaview)
    mean_het = mean(hetview)
    var_het = var(hetview)
    mean_hom = mean(homview)
    var_hom = var(homview)
    #
    mean_a = mean(aview)
    var_a = var(aview)
    mean_d = mean(dview)
    var_d = var(dview)
    submodel.output_stats[i,:] = [maf submodel.k mean_a var_a mean_d var_d 0 0 0 0 mean_i mean_alpha var_alpha 0 0 mean_het var_het mean_hom var_hom 0 0 0 0]
end

function calc_additive_only_stats!(submodel, testsnp, i)
    X = [ones(submodel.n) testsnp.genotypes[:]]
    invXTX = 0
    try
        invXTX = inv(X' * (submodel.weights.rinverse .* X))
    catch
        submodel.output_stats[i,:] = fill(-9, 6)
        return
    end
    XTY = (submodel.weights.rinverse .* X)' * submodel.phenotype.trait
    betahat = invXTX * XTY
    mse = diag(submodel.yty .- betahat'XTY) ./ (submodel.n-2)
    sigma2 = diag(invXTX)*mse'
    sigma = sqrt.(sigma2)
    z = randn(submodel.rng, Float64, (2,submodel.k))
    beta = (z .* sigma) .+ betahat # 3 x k
    if submodel.printchain
        writedlm(submodel.output_filename .* "." .* submodel.snpset.markerInfo[i,3] .* ".tsv", beta_plaus', '\t')
    end
    maf = sum(Array{Int64,1}(testsnp.genotypes[:]))/(2*submodel.n)
    aview = view(beta, 2, :)
    mean_a = mean(aview)
    var_a = var(aview)
    submodel.output_stats[i,:] = [maf submodel.k mean_a var_a 0 0]
end

