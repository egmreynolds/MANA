# Functions of reading phenotype, weights files and adding them to the model
#Read adjusted phenotypes
function read_yk(ykfile::String)
    yk = readdlm(ykfile, '\t', Float64, '\n'; header = false)
    return Phenotypes(map(x->string(Int64(x)), yk[:,1]), yk[:,2:end])
end
# LOSO adjusted
function setup_phenotypes(model, prefix, weightsfilename)
    model.phenotypesprefix = prefix
    model.phenotypesfilename =  string(model.phenotypesprefix, ".layer1.n1.tsv")
    model.phenotypes = read_yk(model.phenotypesfilename)

    model.weightsfilename = weightsfilename
    model.weights = read_weights(model.weightsfilename)

    model.phenoOrder, model.weightOrder = get_orders(model.phenotypes.IID, model.weights.IID)
    get_phenotype_intersect!(model.phenotypes, model.phenoOrder)
    get_weight_intersect!(model.weights, model.weightOrder)
end
#File with IID, Rinv, same order as phenotype file
function read_weights(weightsfile::String)
    w  = readdlm(weightsfile, ' ', Float64, '\n'; header = false)
    return Weights(map(x->string(Int64(x)), w[:,1]), w[:,2])
end
# Full MA adjusted (used full phenotype filename, not prefix)
function setup_phenotypes_full(model, phenotypesfilename, weightsfilename)
    model.phenotypesfilename =  phenotypesfilename
    model.phenotypes = read_yk(model.phenotypesfilename)

    model.weightsfilename = weightsfilename
    model.weights = read_weights(model.weightsfilename)

    model.phenoOrder, model.weightOrder = get_orders(model.phenotypes.IID, model.weights.IID)
    get_phenotype_intersect!(model.phenotypes, model.phenoOrder)
    get_weight_intersect!(model.weights, model.weightOrder)
end

function read_covariates(covarfile)
    covariates = readdlm(covarfile, ' ', Float64, '\n'; header = false)
    return Covariates(map(x->string(Int64(x)), covariates[:,1]), covariates[:,2:end])
end