function [psd, f] = getPSD(spectra, state)
    f = spectra.state(state).f;
    psdMat = spectra.state(state).psd;
    psd = zeros(size(psdMat, 1), size(psdMat, 2));
    for i = 1:size(psdMat, 1)
        psd(i,:) = diag(squeeze(psdMat(i,:,:)));
    end
end