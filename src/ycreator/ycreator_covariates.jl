#!/bin/julia
# functions to adjust Yk phenotypes by a pre-generated set of covariates.

function covariate_setup(phenofile, covarfile)
    phenotypes = read_yk(phenofile)
    covariates = read_covariates(covarfile)
    phenoOrder, covarOrder = get_orders(phenotypes.IID, covariates.IID)
    get_phenotype_intersect!(phenotypes, phenoOrder)
    get_covariate_intersect!(covariates, covarOrder)
    return Phenotypes(phenotypes.IID, calculate_residual_phenotype(phenotypes.trait, covariates.X))
end

function calculate_residual_phenotype(Y, X)
    invXX = inv(X'X)
    XinvXX = X * invXX
    XY = X'Y
    XB = XinvXX * XY
    resY = Y .- XB
    return resY    
end
