# Functions for reading and using plink binary files (.bed, .bim, .fam)
function read_bed_file(bedfile)
    bf = BEDFile(BEDFiles.datadir(bedfile))
    counts(bf, dims=1)
    counts(bf, dims=2)
    return bf
end

function read_bim_file(bimfile)
    bim = readdlm(bimfile, '\t', String, '\n')
    return Array{String,2}(view(bim,:,[1,4,2,5,6]))  ## order it Chr, BP, snp, A1, A2
end

function read_fam_file(famfile)
    fam = readdlm(famfile, '\t', String, '\n')
    try
        return fam[:,2]
    catch
        println("Fam file is not <tab> delimited... trying <space> delimiter (error if <space> also fails)")
        fam = readdlm(famfile, ' ', String, '\n')
    end
    return fam[:,2]
end
# read binary file set
function read_binary_genotypes(fileprefix)
    info = read_bim_file(fileprefix .* ".bim")
    iid = read_fam_file(fileprefix .* ".fam")
    bed = read_bed_file(fileprefix .* ".bed") # remove @time
    return BinaryMarkers(info, iid, bed)
end
#
function add_snps(model, fileprefix)
    model.snpsfile = fileprefix
    model.snpsOld = read_binary_genotypes(fileprefix)
end
