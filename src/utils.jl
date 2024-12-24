function dict_not_subset(
    dict1::Dict{String, Dict{String, Vector{String}}},
    dict2::Dict{String, Dict{String, Vector{String}}})::Bool
    for (key, inner_dict1) in dict1
        if !(key in dict2)
            return true
        end

        inner_dict2 = dict2[key]
        for (inner_key, values1) in inner_dict1
            if !(inner_key in inner_dict2)
                return true
            end
            values2 = inner_dict2[inner_key]
                if !(values1 ⊆ values2)
                    return true
                end
        end
    end
    return false
end