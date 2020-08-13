function [vpath, Gamma] = processVpathAMICA(v)
    rej_index = any(sum(v, 1)==0, 1);
    valid_index = logical(1 - rej_index);
  
    v(:,valid_index) = 10.^v(:, valid_index);
    
    Gamma = v';
    
    [~, vpath] = max(v, [], 1);
    vpath(:,rej_index) = 0;
    vpath = vpath';
end