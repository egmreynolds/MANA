# Functions assisting those in run_base.jl
# get Marker object from BinaryMarker object, at start and finish indices (m: markers, i:indivs)
function get_bed_values(markers::BinMarkers, m1, m2)
    genotypes = bedvals[view(markers.binarygenotypes, :, m1:m2)]
    flip_alleles!(genotypes)
    info = markers.markerInfo[m1:m2,:]
    return Markers(info, markers.IID, Array{Float64,2}(genotypes))
end

function calc_betas(marker, y::Array{Float64,2}, n::Int64, rinv)
    m = Array{Int64,1}(marker .+ 1)
    X = sparse(1:n, m, 1.0)
    rinvX = rinv .* X
    cellcounts = sum(rinvX, dims=1)
    if cells_zero(cellcounts')
        return zeros(3,1), zeros(3,1), true
    end
    cellsums = rinvX'y
    return cellsums, cellcounts', false
end

function select_snp(submodel, i)
    return Markers(submodel.snpset.markerInfo[i:i,:], submodel.snpset.IID, submodel.snpset.genotypes[:,i:i])
end

function get_1loci_k()
    return Array{Float64,2}([1.0 0.0 0.0
    -0.5 0.0 0.5
    -0.5 1.0 -0.5])
end

function cells_zero(cells)
    return (size(cells,1) < 3 || cells[1] <= 5 || cells[2] <= 5 || cells[3] <= 5)
end

function next_phenotype_interval!(model)
    model.layerlimit += 0.5*model.losoInterval
    println("New Limit: ", model.layerlimit)
    if model.layerlimit < model.maxlayerlimit
        if model.inlayer1
            model.nlayer2 += 1
            model.phenotypesfilename = string(model.phenotypesprefix, ".layer2.n", model.nlayer2, ".tsv")
            model.phenotypes = read_yk(model.phenotypesfilename)
        else
            model.nlayer1 += 1
            model.phenotypesfilename = string(model.phenotypesprefix, ".layer1.n", model.nlayer1, ".tsv")
            model.phenotypes = read_yk(model.phenotypesfilename)
        end
        model.inlayer1 = !model.inlayer1
        get_phenotype_intersect!(model.phenotypes, model.phenoOrder)
        get_phenotype_intersect!(model.phenotypes, model.phenoOrder2)
    end
    println("Switched to new YLOSO interval: ", model.phenotypesfilename)
end
# Most variants seems to have 2 as major allele, 0 as minor, here we flip them.
function flip_alleles!(X)
    for t in eachindex(X)
        if X[t] == 0
            X[t] = 2
        elseif X[t] == 2
            X[t] = 0
        end
    end
end