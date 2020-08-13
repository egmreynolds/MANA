# Functions for adjusting phenotypes (ycreator)
function create_yk_loso(ycreator)
    start = 0
    for i in 1:ycreator.nlayer1
        ycreator.alphas_adj = deepcopy(ycreator.alphas)
        finish = start + ycreator.interval
        restrict_alphas_loso!(ycreator, start, finish)
        Yk = Phenotypes(ycreator.phenotypes.IID, calc_yk(ycreator))
        writedlm(string(ycreator.outputprefix, ".", Int64(ycreator.interval/1000000), "Mb.layer1.n", Int64(i), ".tsv"), hcat(Yk.IID, Yk.trait))
        start = finish
    end
    start = ycreator.interval ./ 2.0
    for i in 1:ycreator.nlayer2
        ycreator.alphas_adj = deepcopy(ycreator.alphas)
        finish = start + ycreator.interval
        restrict_alphas_loso!(ycreator, start, finish)
        Yk = Phenotypes(ycreator.phenotypes.IID, calc_yk(ycreator))
        writedlm(string(ycreator.outputprefix, ".", Int64(ycreator.interval/1000000), "Mb.layer2.n", Int64(i), ".tsv"), hcat(Yk.IID, Yk.trait))
        start = finish
    end
end

function restrict_alphas_loso!(ycreator, start, finish)
    chridx = findall(x -> x == ycreator.chr, ycreator.genselmarkers.SNPID[:,1])
    bpidx = findall(x -> (parse(Float64,x) >= start && parse(Float64,x) <= finish), ycreator.genselmarkers.SNPID[:,2])
    intersectIndex = intersect(chridx, bpidx)
    ycreator.alphas_adj.alphas[intersectIndex, :] .= 0.0
end

function calc_yk(ycreator)  # y:nx1, m:nxm, k:mxk
    ma = Array{Float64,2}(ycreator.genselmarkers.markers) * Array{Float64,2}(ycreator.alphas_adj.alphas)
    return ycreator.phenotypes.trait .- ma
end

## Additional functions to general Ma intermediate step. Modifiers on above functions..
function create_Ma_loso(ycreator)
    start = 0
    for i in 1:ycreator.nlayer1
        ycreator.alphas_adj = deepcopy(ycreator.alphas)
        finish = start + ycreator.interval
        restrict_alphas_loso!(ycreator, start, finish)
        Ma = calculate_Ma(ycreator)
        writedlm(string(ycreator.outputprefix, ".", Int64(ycreator.interval/1000000), "Mb.layer1.n", Int64(i), ".tsv"), hcat(ycreator.phenotypes.IID, Ma))
        start = finish
    end
    start = ycreator.interval ./ 2.0
    for i in 1:ycreator.nlayer2
        ycreator.alphas_adj = deepcopy(ycreator.alphas)
        finish = start + ycreator.interval
        restrict_alphas_loso!(ycreator, start, finish)
        Ma = calculate_Ma(ycreator)
        writedlm(string(ycreator.outputprefix, ".", Int64(ycreator.interval/1000000), "Mb.layer2.n", Int64(i), ".tsv"), hcat(ycreator.phenotypes.IID, Ma))
        start = finish
    end
end
function calculate_Ma(ycreator)
    ma = Array{Float64,2}(ycreator.genselmarkers.markers) * Array{Float64,2}(ycreator.alphas_adj.alphas)
    return ma
end

