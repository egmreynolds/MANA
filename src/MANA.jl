# Packge Main for MANA GWAS package
using Distributed, DelimitedFiles, Statistics, BGZFStreams, LinearAlgebra, Random, SparseArrays, SharedArrays, BEDFiles, Distributions
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
#
export build_gwas
export add_snps
export setup_phenotypes
export initiate_model
export run_model
export concatenate_summaries
#
export make_alphas
export filter_alphas
export print_alphas
#
export build_ycreator
export add_phenotypes
export add_alphas
export add_markers_gensel
export sort_ycreator
export initiate_ycreator
export create_yk_loso
export create_yk_full
export run_2locus_model
export setup_phenotypes_full
export run_pedigree_model
export run_spliceQTL 
export run_pedigree_model_multirecord
