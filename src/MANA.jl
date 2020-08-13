module MANA
# Packge Main for MANA GWAS package
using Distributed, DelimitedFiles, Statistics, LinearAlgebra, Random, SparseArrays, BEDFiles, Distributions
using JuliaDB, IndexedTables, OnlineStats, JWAS, DataFrames
#
include("types.jl")
include("read/read.jl")
include("model/model.jl")
include("mutate/mutate.jl")
include("run/run.jl")
include("summarise/summarise.jl")
include("ycreator/ycreator.jl")
include("pedigree/pedigree.jl")

# Functions for basic model
export build_gwas, add_snps, setup_phenotypes, initiate_model, run_model, concatenate_summaries,
# Functions for adjusting marker effects
export make_alphas, filter_alphas, print_alphas
# Functions for adjusting phenotypes
export build_ycreator, add_phenotypes, add_alphas, add_markers_gensel, sort_ycreator, initiate_ycreator, create_yk_loso
# Functions for pedigree models
export run_pedigree_model, run_spliceQTL, run_pedigree_model_multirecord

end
