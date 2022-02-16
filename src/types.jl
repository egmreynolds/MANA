mutable struct Model
    snpsOld
    snpsfile
    phenotypesprefix
    phenotypesfilename
    phenotypes
    m
    maxbp
    nlayer1
    nlayer2
    losoInterval
    layerlimit
    inlayer1
    setsize
    notcomplete
    phenoOrder
    phenoOrder2
    genoOrder
    snps
    n
    k
    outputfilename
    printchain
    submodelcount
    futures
    weightsfilename
    weights
    weightOrder
    maxlayerlimit
    type
end

mutable struct SubModel
    m
    n
    k
    snpset
    phenotype
    phenotypesfilename
    K
    rng
    yty
    output_stats
    output_filename
    printchain
    id
    weights
    type
end

mutable struct BinaryMarkers{M<:Array{String,2}, I<:Array{String,1}, B<:BEDFile}
    markerInfo::M
    IID::I
    binarygenotypes::B
end

mutable struct Phenotypes{I<:Array{String,1}, T}
    IID::I
    trait::T # n x k
end

mutable struct Weights{I<:Array{String,1}, W<:Array{Float64,1}}
    IID::I
    rinverse::W
end


mutable struct BinMarkers{M<:Array{String,2}, I<:Array{String,1}, B<:Array{UInt8,2}}
    markerInfo::M
    IID::I
    binarygenotypes::B
end

mutable struct Markers{M<:Array{String,2}, I<:Array{String,1}, G<:Array{Float64,2}}
    markerInfo::M # m x 5 columns
    IID::I # n x 1
    genotypes::G # n x m
end

mutable struct Alphas
    markerNames
    alphas
    m
    k
end

mutable struct YCreator
    phenotypes
    alphas
    genselmarkers
    n
    m
    k
    ynew
    nlayer1
    nlayer2
    alphas_adj
    outputprefix
    interval
    chr
    chrlength
    region
    code
end

mutable struct GenselMarkers
    IID
    SNPID
    markers
end

mutable struct id_dataset
        iid
        header
        data
end

mutable struct Covariates{I<:Array{String,1}, x<:Array{Float64,2}}
    IID::I
    X::x
end
