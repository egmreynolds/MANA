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
export build_gwas
export add_snps
export setup_phenotypes
export initiate_model
export run_model
export concatenate_summaries
# Functions for adjusting marker effects
export make_alphas
export filter_alphas
export print_alphas
# Functions for adjusting phenotypes
export build_ycreator
export add_phenotypes
export add_alphas
export add_markers_gensel
export sort_ycreator
export initiate_ycreator
export create_yk_loso
# Functions for pedigree models
export run_pedigree_model
export run_spliceQTL
export run_pedigree_model_multirecord

end
