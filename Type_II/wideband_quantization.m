function [P,W_max] = wideband_quantization(W)
W = W(:);
W = abs(W);
W_max = max(W);
W_order = W/W_max;
codebook = [1, 1/sqrt(2), 1/sqrt(4), 1/sqrt(8), 1/sqrt(16), 1/sqrt(32), 1/sqrt(64), 0];
P = quantize(W_order, codebook);
P = diag(P);
end

function ret = quantize(samples, codes)
    codes = codes(:);
    codes = sort(codes)';
    midpoints = (codes(2 : end) + codes(1 : end-1)) / 2;
    diffs = codes(2 : end) - codes(1 : end-1);
    ret = sum((samples(:) > midpoints) .* diffs, 2) + codes(1);
    ret = reshape(ret, size(samples));
end
