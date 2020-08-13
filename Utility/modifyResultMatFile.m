
fields = {'GammaInit', 'Xi'};

for K = 2:10
    K
    cd(strcat('/home/ting/Documents/eeg_hmm/HMM_results/results_K', num2str(K)));
    datafile_list = dir('*.mat');
    
    for i = 1:length(datafile_list)
        filename = datafile_list(i).name
        try
            data = load(filename);
        catch error
            disp(strcat('error loading', {' '}, filename));
            continue
        end
        
        for field = fields
            if isfield(data, field{1})
                data = rmfield(data, field{1});
            end
        end
        save(filename, '-struct', 'data');
        clear 'data'
    end
end



