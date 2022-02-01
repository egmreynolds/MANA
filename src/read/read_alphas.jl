# Functions for reading, using, analysing gensel outputs (alphas)
# alphas come in binary format (m=Int32, m*{Int32, Float32})
function read_bin_alphas(alphas_filename, m)
    stream = open(alphas_filename, "r")
    alphas = Array{Float64}(1:m)
    while !eof(stream)
        mx = read(stream, Int32)
        alphas_col = zeros(Float64, m)
        for i in 1:mx
            index = read(stream, Int32) + 1
            alpha = read(stream, Float32)
            alphas_col[index] = alpha
        end
        alphas = hcat(alphas, alphas_col)
    end
    return alphas
end
# mrkRes files
function read_alphanames(alphanames_filename)
    alpharesults = readdlm(alphanames_filename, ' ', String, '\n', header = true)
    return alpharesults[1][:,1]
end
# keep every 'freq'th sample between 'first' and 'last' index.
# Note every 'freq'th means if modulo freq == 0. NOT from first.
function filter_alphas(alphas, freq, first, final)
    filterindex = Int64[]
    for i in first:final
        if (i % freq) == 0
            push!(filterindex, i)
        end
    end
    alphas.alphas = alphas.alphas[:,filterindex]
    alphas.k = size(alphas.alphas, 2)
end

function print_alphas(alphas, output_filename)
    output_alphas = hcat(alphas.markerNames, alphas.alphas)
    writedlm(output_filename, output_alphas)
end
# Make Alphas type
function make_alphas(alphas, alphanames_filename, alphas_filename)
    alphas.markerNames = read_alphanames(alphanames_filename)
    alphas.m = size(alphas.markerNames, 1)
    alphas.alphas = read_bin_alphas(alphas_filename, alphas.m)
    alphas.k = size(alphas.alphas, 2)
end
# geweke diagnostics of MCMC chains
function get_mcmc_diagnostics(alphas, outputprefix)
    geweke_pvalues = zeros(Float64, alphas.m)
    try 
        geweke_pvalues = gewekediag(Chains(alphas.alphas[:,2:end]'))[1][:,3]
    catch        
        for i in 1:alphas.m
            geweke_pvalues[i] = gewekediag(alphas.alphas[i,2:end])[2]
        end        
    end
    geweke_by_marker = hcat(alphas.markerNames, geweke_pvalues)
    writedlm(outputprefix .* ".geweke.stats", geweke_by_marker)

    geweke_pass = sum(geweke_pvalues .> 0.05)
    writedlm(outputprefix .* ".geweke.summary", hcat(alphas.m, geweke_pass, geweke_pass/alphas.m))
end
# Initiate an Alphas object
function build_alphas()
    return Alphas(0,0,0,0)
end
