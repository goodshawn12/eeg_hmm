function transProb = getTransProbMatrix(vpath, K)
    transProb = zeros(K);
    
    invalid = find(~vpath);
    if ~isempty(invalid)
        % The first index of invalid cannot be -1
        if invalid(1) == 1  
            invalid(1) = [];
        end    
        invalid = union(invalid - 1, invalid);
    end    
    % The last index is always invalid
    invalid = [invalid; length(vpath)];
    
    valid = setdiff(1:length(vpath), invalid);
    
    count_arr = zeros(K, 1);    
    for prev = valid
        oldState = vpath(prev);
        newState = vpath(prev+1);
        transProb(oldState, newState) = 1 + transProb(oldState, newState);
        count_arr(oldState) = count_arr(oldState) + 1;
    end
    transProb = transProb ./ count_arr;    
end