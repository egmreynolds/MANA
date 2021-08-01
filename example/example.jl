using MANA

function main()
	dir = "~/MANA/example/" # May need to be edited to directory of interest
	snpfilename = dir .* "plink.example"
        phenotypesprefix = dir .* "example1_phenotypes/phen.example"
        outputfilename = dir .* "out.example"
        printchain = false
        weightsfilename= dir .* "weights.example.txt"
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

main()
