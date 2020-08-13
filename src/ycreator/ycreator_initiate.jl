# Functions used to initiate ycreator type, for adjusting phenotypes with gensel outputs
function build_ycreator()
    return YCreator(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
end
#
function add_markers_gensel(ycreator, markersfileprefix)
    info = read_bim_file(markersfileprefix .* ".bim")
    iid = read_fam_file(markersfileprefix .* ".fam")
    bed = read_bed_file(markersfileprefix .* ".bed")
    genselvalues = (bedvals[bed] .- 1).*-1  # convert to -1, 0, 1 # *-1 as plinks ref/alt coding is differnt for bed files vs raw files. gensel alphas are calculated wrt to raw files.
    ycreator.genselmarkers = GenselMarkers(iid, info, genselvalues)
end
#
function add_alphas(ycreator, alphasfilename)
    read_filtered_alphas(ycreator, alphasfilename)
end
# rsplit acts in reverse, and limit means it just splits on the 1st occurence. (This prevents removal of missingIDs (CHRN_POS) format)
function read_filtered_alphas(ycreator, alphasfilename)
    alphas = readdlm(alphasfilename, '\t', Any, '\n'; header = false)
    ycreator.alphas = Alphas(alphas[:,1], alphas[:,2:end], 0 ,0)
    ycreator.alphas.m = size(ycreator.alphas.markerNames, 1)
    ycreator.alphas.k = size(ycreator.alphas.alphas, 2)
    ycreator.alphas.markerNames = map(x -> String(x[1]), rsplit.(ycreator.alphas.markerNames, '_'; limit = 2))
end

function add_phenotypes(ycreator, phenotypesfilename)
    ycreator.phenotypes = read_y(phenotypesfilename)

end

function read_y(phenotypesfilename)
    y = readdlm(phenotypesfilename, ' ', Float64, '\n'; header = true)
    return Phenotypes(map(x->string(Int64(x)), y[1][:,1]), y[1][:,2:2])
end
# marker and alpha orders should be the same
function sort_ycreator(ycreator)
    markerOrderM = collect(1:ycreator.alphas.m)
    markerOrderN, phenoOrderN = get_orders(ycreator.genselmarkers.IID, ycreator.phenotypes.IID)
    get_genselmarker_intersect!(ycreator.genselmarkers, markerOrderN, markerOrderM)
    get_phenotype_intersect!(ycreator.phenotypes, phenoOrderN)
end

function get_genselmarker_intersect!(genselmarkers, new_orderN, new_orderM)
    genselmarkers.IID = genselmarkers.IID[new_orderN]
    genselmarkers.SNPID = genselmarkers.SNPID[new_orderM, :]
    genselmarkers.markers = genselmarkers.markers[new_orderN, new_orderM]
end

function get_alpha_intersect!(alphas, new_order)
    alphas.markerNames = alphas.markerNames[new_order]
    alphas.alphas = alphas.alphas[new_order, :]
    alphas.m = size(alphas.markerNames, 1)
end

function initiate_ycreator(ycreator, interval, chr, outputprefix, region, code)
    ycreator.n = size(ycreator.phenotypes.IID, 1)
    ycreator.m = size(ycreator.genselmarkers.SNPID, 1)
    ycreator.k = size(ycreator.alphas.k)
    ycreator.outputprefix = outputprefix
    ycreator.interval = interval
    ycreator.region = region
    ycreator.chr = chr
    ycreator.code = code
    if ycreator.code == "LOCO"
        ycreator.chrlength = 0
    else
        ycreator.chrlength = get_chr_length(ycreator.genselmarkers.SNPID, ycreator.chr)
    end
    ycreator.nlayer1 = ceil(ycreator.chrlength / ycreator.interval)
    ycreator.nlayer2 = ycreator.nlayer1 - 1
end

function get_chr_length(markerinfo, chr)
    chrindex = findall(x -> x == chr, markerinfo[:,1])
    chrview = view(markerinfo, chrindex, 2)
    chrviewfloat = map(x -> parse(Float64, x), chrview)
    chrlength = maximum(chrviewfloat)
    return chrlength
end

function fill_missing_genotypes!(ycreator)
    if any(ismissing.(ycreator.genselmarkers.markers))
        ycreator.genselmarkers.markers = Array{Union{Missing, Float64},2}(ycreator.genselmarkers.markers)
        for i in 1:ycreator.m
            l = length(ycreator.genselmarkers.markers[ismissing.(ycreator.genselmarkers.markers[:,i]),i])
            if l != 0
                ycreator.genselmarkers.markers[ismissing.(ycreator.genselmarkers.markers[:,i]),i] = fill(mean(skipmissing(ycreator.genselmarkers.markers[:,i])),l)
            end
        end
    end
end
