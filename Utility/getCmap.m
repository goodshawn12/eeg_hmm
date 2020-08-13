function [cmap, state_description] = getCmap(K)
    if K == 5
        cmap = [1,0.1,0.1; 1,0.8,0.2; 0.8,0.7,0.55; 0.1,1,0.1; 0.1,0.1,1];
        state_description = {'Drowsy', 'Middle1', 'Middle2', 'Middle3', 'Alert'};
    elseif K == 4
        cmap = [1,0.1,0.1; 1,0.8,0.2; 0.1,1,0.1; 0.1,0.1,1];
        state_description = {'Drowsy', 'Middle1', 'Middle2', 'Alert'};
    elseif K == 3
        cmap = [1,0.1,0.1; 0.1,1,0.1; 0.1,0.1,1];
        state_description = {'Drowsy', 'Middle1', 'Alert'};
    elseif K == 2
        cmap = [1,0.1,0.1; 0.1,0.1,1];
        state_description = {'Drowsy', 'Alert'};
    else
        cmap = parula;
    end
end