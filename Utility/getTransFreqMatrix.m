% state-state transition count divided by number of transitions (including identity transition)
function transFreq = getTransFreqMatrix(vpath, K)
    transFreq = zeros(K);
    
    invalid = find(~vpath);
    % The first index of invalid cannot be -1
    if invalid(1) == 1  
        invalid(1) = [];
    end    
    invalid = union(invalid - 1, invalid);
    % The last index is always invalid
    invalid = [invalid; length(vpath)];
    
    valid = setdiff(1:length(vpath), invalid);
    
    for prev = valid
        oldState = vpath(prev);
        newState = vpath(prev+1);
        transFreq(oldState, newState) = 1 + transFreq(oldState, newState);
    end
    transFreq = transFreq ./ length(valid);
end