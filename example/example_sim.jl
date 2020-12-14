# Large example based on simulated data
# Create LOSO adjusted phenotypes based on GenSel output
# Run GWAS based on these
#
# Download dataset at https://datadryad.org/ 
# Extract:
# 800K genotypes on 80K indivs: gen11.1.bed/bim/fam
# 31K genotypes on 80K indivs: gen12.1.bed/bim/fam
# Simulated phenotype: phen_8.1.txt
#
# These two functions will typically be run in separate chunks as they will typically have different CPU/MEM requirements
#
using MANA
# create adjusted phenotypes
function make_loso_phenotypes()
  chr=26
  interval=10000000
  phenotypesfilename=phen_8.1.txt
  markersfileprefix=gen12.1
  alphasfilename=other_13.txt
  outputprefix=adjusted_phenotype.chr26
  REGION=0
  CODE="LOSO"
  ycreator = build_ycreator()
  add_alphas(ycreator, alphasfilename) # marker effects
  add_phenotypes(ycreator, phenotypesfilename) # phenotype
  add_markers_gensel(ycreator, markersfileprefix) # marker genotypes
  sort_ycreator(ycreator) 
  initiate_ycreator(ycreator, interval, chr, outputprefix, REGION, CODE)  # set up creator object
  create_yk_loso(ycreator) # create adjusted phenotypes based on LOSO and the interval given
end
#
make_loso_phenotype()
#
# run GWAS
function main()
	dir = "/path/to/files/" # May need to be edited to directory of interest
	snpfilename = dir .* "gen11.1"
  phenotypesprefix = dir .* "adjusted_phenotype.chr26"
  outputfilename = dir .* "out.simulation_example"
  printchain = false
  weightsfilename= dir .* "weights.simulation.txt" # file of IDs and weights (default = 1)
  interval = 10000000
  setsize = 1000
  type = "dominance"
  # Build Model
  model = build_gwas()
  add_snps(model, snpfilename)
  setup_phenotypes(model, phenotypesprefix, weightsfilename)
  initiate_model(model, interval, setsize, outputfilename, printchain, type)
  # Run Model
  run_model(model)
  # Summarise
  concatenate_summaries(model)
end
#
main()
