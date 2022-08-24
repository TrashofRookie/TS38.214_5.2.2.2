function [quan,quan_polar] = wideband_quantization(b_max,b_max_polar)
W = [b_max b_max_polar];
W = abs(W);
W_max = max(W);
W = W / W_max;
codebook = [1, (1/2)^(1/4), 1/sqrt(2), (1/8)^(1/4), 1/2, (1/32)^(1/4), 1/sqrt(8), (1/128)^(1/4), 1/4, (1/512)^(1/4), 1/(2 * sqrt(8)), (1/2048)^(1/4), 1/8, (1/8192)^(1/4), 1/sqrt(128), 0];
P = quantize(W, codebook);
quan = P(1);
quan_polar = P(2);
end

function ret = quantize(samples, codes)
    codes = codes(:);
    codes = sort(codes)';
    midpoints = (codes(2 : end) + codes(1 : end-1)) / 2;
    diffs = codes(2 : end) - codes(1 : end-1);
    ret = sum((samples(:) > midpoints) .* diffs, 2) + codes(1);
    ret = reshape(ret, size(samples));
end
