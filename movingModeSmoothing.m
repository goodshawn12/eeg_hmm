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
    for index = 1:size(A, mod(dim,2)+1)
        i = half_window + 1;
        while i < size(A,dim)
            window_start = i - half_window;
            window_end = i + half_window;
            if window_end > size(A, dim)
                window_end = size(A, dim);
            end
            
            if dim == 1
                B(i, index) = quickMode(A(window_start:window_end, index), values);
            else
                B(index, i) = quickMode(A(index, window_start:window_end), values);
            end            
            i = i + 1;
        end
    end
end