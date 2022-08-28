function [P] = narrowband_quantization(W)
W = abs(W);
W_max = max(W);
W_norm = W / W_max;
codebook = [1, 1/sqrt(2), 1/2, 1/sqrt(8), 1/4, 1/sqrt(32), 1/8, 1/sqrt(128)];
P = quantize(W_norm, codebook);
 
end

function ret = quantize(samples, codes)
    codes = codes(:);
    codes = sort(codes)';
    midpoints = (codes(2 : end) + codes(1 : end - 1)) / 2;
    diffs = codes(2 : end) - codes(1 : end-1);
    ret = sum((samples(:) > midpoints) .* diffs, 2) + codes(1);
    ret = reshape(ret, size(samples));
end
