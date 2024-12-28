function dict_not_subset(dict1::Dict, dict2::Dict)::Bool
    return !all(k -> haskey(dict2, k) && sort(dict1[k]) == sort(dict2[k]), keys(dict1))
end