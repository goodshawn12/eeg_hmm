function [merged_struct] = mergeStructs(struct1, struct2)
    if isempty(struct1)
        merged_struct = struct2;
        return
    end
    if isempty(struct2)
        merged_struct = struct1;
        return
    end
    
    merged_struct = struct1;
    f = fieldnames(struct2);
    for i = 1:length(f)
        merged_struct.(f{i}) = struct2.(f{i});
    end
end