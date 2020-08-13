# Functions focusing on finding relative orderings between datasets and reordering those datasets accordingly
function get_orders(order1::Array{String,1}, order2::Array{String,1})
    common_ids = intersect(order1, order2)
    new_order1 = get_id_positions(common_ids, order1)
    new_order2 = get_id_positions(common_ids, order2)
    return new_order1, new_order2
end
#
function get_id_positions(ids::Array{String,1}, old_order::Array{String,1})
    new_order = zeros(Int64, size(ids,1))
    for i in 1:size(ids,1)
        for j in 1:size(old_order,1)
            if ids[i] == old_order[j]
                new_order[i] = j
                break
            end
        end
    end
    return new_order
end
#
function get_phenotype_intersect!(phenotypes::Phenotypes, new_order::Array{Int64,1})
    phenotypes.trait = phenotypes.trait[new_order,:]
    phenotypes.IID = phenotypes.IID[new_order]
end
# Orders snps and also converts to a BinMarker object from BinaryMarkers
function get_genotype_intersect_bin(markers::BinaryMarkers, new_order::Array{Int64,1})
    return BinMarkers(markers.markerInfo, markers.IID[new_order], markers.binarygenotypes[new_order, :])
end

function get_weight_intersect!(weights::Weights, new_order::Array{Int64,1})
    weights.rinverse = weights.rinverse[new_order]
    weights.IID = weights.IID[new_order]
end
