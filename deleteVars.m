function result_table = deleteVars(table, vars)
    % Change single string to cell
    if ischar(vars)
        vars = {vars};
    end
    % Throw error if not cell or table
    if ~iscell(vars)
        error('Only cell array of strings or single string')
    end
    if ~istable(table)
        error('First argument must be a table')
    end
    % Throw error if cell has more dimensions
    if length(size(vars)) > 2 || (size(vars, 1) ~= 1 && size(vars, 2) ~= 1)
        error('Only accept a 1D cell array')
    end
    
    for i = 1:length(vars)
        if ismember(vars{i}, table.Properties.VariableNames)
            table = removevars(table, vars{i});
        end
    end
    
    result_table = table;
end