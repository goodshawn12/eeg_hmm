function B = movingModeSmoothing(A, window_size, dim, values)
    if window_size < 3
        B = A;
        return
    end
    if length(size(A)) > 2 || dim > 2
        error('Only support 2-d matrix')
    end
    if mod(window_size, 2) == 0
        warning('Even window length will be treated as window_size - 1')
    end
    B = A;
    half_window = floor(window_size/2);
    
    for index = 1:size(A, 3-dim)
        j = half_window + 1;
        slice = zeros(1, size(A, dim));
        while j < size(A,dim)
            window_start = j - half_window;
            window_end = window_start + window_size - 1;
            if window_end > size(A, dim)
                window_end = size(A, dim);
            end
            
            if dim == 1
                B(j, index) = quickMode(A(window_start:window_end, index), values);
            else
                B(index, j) = quickMode(A(index, window_start:window_end), values);
            end            
            j = j + 1;
        end
    end
    delete(gcp('nocreate')); % shut down any current pool
end