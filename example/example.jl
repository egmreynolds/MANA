using MANA

function main()
	dir = "/data/seq/edrey0/MANA/example/"
	snpfilename = dir .* "plink.example"
        phenotypesprefix = dir .* "phen.example"
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
