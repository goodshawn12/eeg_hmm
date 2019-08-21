function B = verticalSmoothing(A, window_size)
    B = A;
    for column = 1:size(A, 2)
        i = floor(window_size/2);
        while i < size(A,1)
            window_start = i - floor(window_size/2);
            window_end = i + ceil(window_size/2) - 1;
            if window_start < 1
                window_start = 1;
            end
            if window_end > size(A, 1)
                window_end = size(A, 1);
            end
            B(window_start:window_end, column) = mode(A(window_start:window_end,column));
            i = i + window_size;
        end
    end
end