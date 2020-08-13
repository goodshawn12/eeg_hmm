% For n-by-m data, data is split into n/entries-by-m matrices, and
% the output is npairings-by-rept, where rept is n/nentries (m if dim=2) where 
% npairings is the number of non-identity pairings of 1:m (1:rept if dim=2), i.e., (m-1)+(m-2)+...+1
function [p, h, stats] = mult_pair_signrank(data, n_entries, dim)
    if mod(size(data,1), n_entries) ~= 0
        error('Data along the selected dimension cannot be split evenly by rept times')
    end
%     if size(data, 3-dim) < 2
%         eroor('Need at least two groups to do paired sign rank test')
%     end
    
    if dim == 1
        rept = size(data, 1) / n_entries;
        n_groups = size(data, 2);
    elseif dim == 2
        rept = size(data, 2);
        n_groups = size(data, 1) / n_entries;
    end
    [x_coor, y_coor] = meshgrid(1:n_groups, 1:n_groups);
    x_coor = tril(x_coor, -1);
    y_coor = tril(y_coor, -1);
    group_pairs = [x_coor(find(x_coor)), y_coor(find(y_coor))];   
    n_pairings = size(group_pairs,1);
    p = zeros(n_pairings, rept);
    h = zeros(n_pairings, rept);
    stats = cell(n_pairings, rept);
    
    for rept = 1:rept
        for pair_index = 1:n_pairings
            pair = group_pairs(pair_index, :);
            group1 = pair(1); group2 = pair(2);
            if dim == 1
                row_start = 1+(rept-1)*n_entries;
                x = data(row_start:row_start+n_entries-1, group1);
                y = data(row_start:row_start+n_entries-1, group2);
                [p_val,h_val,stats_val] = signrank(x, y);
            elseif dim == 2
                group1_row_start = 1+(group1-1)*n_entries;
                group2_row_start = 1+(group2-1)*n_entries;
                x = data(group1_row_start:group1_row_start+n_entries-1, rept);
                y = data(group2_row_start:group2_row_start+n_entries-1, rept);
                [p_val,h_val,stats_val] = signrank(x, y);
            end
            p(pair_index, rept) = p_val;
            h(pair_index, rept) = h_val;
            stats{pair_index, rept} = stats_val;
        end
    end
    
end
