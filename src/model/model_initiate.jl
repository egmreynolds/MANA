#Functions for building, initiating models and submodels.
function build_gwas()
    return Model(0,0, 0,0, 0, 0, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0)
end

function initiate_model(model, interval, setsize, outputfilename, printchain, type)
    model.nlayer1 = 1
    model.nlayer2 = 0
    model.losoInterval = Float64(interval)
    model.layerlimit = 0.75 * model.losoInterval
    model.inlayer1 = true
    model.setsize = setsize
    model.notcomplete = true
    model.outputfilename = outputfilename
    model.printchain = printchain
    model.submodelcount = 0
    model.futures = Future[]
    #
    try
        model.m = size(model.snpsOld.binarygenotypes, 2)
        model.maxbp = parse(Float64, model.snpsOld.markerInfo[model.m,2])
    catch
        println("Model already initiated. Please build model again.")
    end
    #
    model.phenoOrder2, model.genoOrder = get_orders(model.phenotypes.IID, model.snpsOld.IID) # keep phenoOrder and phenoOrder2 separate as will have to use both when reading new adjusted phenotypes...
    get_phenotype_intersect!(model.phenotypes, model.phenoOrder2)
    model.weightOrder = model.phenoOrder2 # since pheno and weights have already been sorted/ordered. Weights won't need to be ordered again
    get_weight_intersect!(model.weights, model.weightOrder)
    model.snps = get_genotype_intersect_bin(model.snpsOld, model.genoOrder)
    model.snpsOld = nothing
    model.n = size(model.phenotypes.IID,1)
    model.k = size(model.phenotypes.trait,2)
    model.maxlayerlimit = get_maxlayerlimit(model)
    model.type = type
end

function inititate_submodel!(submodel)
    submodel.rng = MersenneTwister()
    submodel.K = get_1loci_k()    
    submodel.output_stats = zeros(Float64,submodel.m, 23)
    if submodel.type == "additive"
        submodel.output_stats = zeros(Float64, submodel.m, 6)
        submodel.K = 0
    end    
    submodel.yty = (submodel.phenotype.trait .* submodel.weights.rinverse)'submodel.phenotype.trait #(y .* rinv)'y' where rinv is nx1 vector - faster than y'Rinv*y, where Rinv is a diag matrix. Same allocations
end

# read directory of yloso phenotype files and find largest, calculate max layer limit
function get_maxlayerlimit(model)
    phenotypedir = dirname(model.phenotypesprefix)
    phenotypefiles = readdir(phenotypedir)
    maxn = 1
    maxlayer = 1
    phenotypefiles = phenotypefiles[occursin.(basename(model.phenotypesprefix), phenotypefiles)]
    for file in phenotypefiles
        n,layer = get_n_and_layer(file)
        if n >  maxn
            maxn = n
            maxlayer = layer
        elseif n == maxn
            if layer > maxlayer
                maxlayer = layer
                maxn = n
            end
        end
    end
    maxlayerlimit = maxn * 10000000
    if maxlayer == 2
        maxlayerlimit = maxlayerlimit + 5000000
    end
    return maxlayerlimit
end

function get_n_and_layer(x)
    sub = x[1:findlast(isequal('.'),x)-1]
    n = parse(Int64, split(sub[(findlast(isequal('.'),sub)+1):end], "n")[2])
    sub2 = sub[1:findlast(isequal('.'),sub)-1]
    layer = parse(Int64, split(sub2[(findlast(isequal('.'),sub2)+1):end], "layer")[2])
    return n, layer
end
