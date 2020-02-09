function [p] = welchttest(X1, X2)
    % Welch's t-test two-tailed
    
    if length(size(X1)) > 2 || length(size(X1)) > 2 || ~any(size(X1)==1) || ~any(size(X2)==1)
        error('Input should be 1-d arrays')
    end
    
    Z = mean(X1) - mean(X2);
    s1 = var(X1)/length(X1);
    s2 = var(X2)/length(X2);
    s = sqrt(s1 + s2);
    t = Z/s;
    
    dof1 = length(X1) - 1;
    dof2 = length(X2) -1;
    dof = (s1 + s2)^2 / ((s1^2)/dof1 + (s2^2)/dof2);
    
    if t < 0
        p = tcdf(t, dof) * 2;
    else
        p = (1 - tcdf(t, dof)) * 2;
    end
end