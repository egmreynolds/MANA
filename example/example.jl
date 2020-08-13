function main()
	snpfilename = "plink.example"
        phenotypesprefix = "phen.example"
        outputfilename = "out.example"
        printchain = false
        weightsfilename= "weights.example.txt"
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
