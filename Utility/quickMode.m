function mode = quickMode(data, values)
%QUICKMODE Summary of this function goes here
%   Detailed explanation goes here
    count = zeros(size(values));
    for i = 1:length(values)
        count(i) = sum(data == values(i));
    end
    [~, index] = max(count);
    mode = values(index);
end

