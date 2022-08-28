function [p_WB] = wideband_quantization(w2_max,w2_max_polar, numberofBeams)
W = [w2_max w2_max_polar];
W = abs(W);
W_max = max(W);
W = W / W_max;
codebook = [1, (1/2)^(1/4), 1/sqrt(2), (1/8)^(1/4), 1/2, (1/32)^(1/4), 1/sqrt(8), (1/128)^(1/4), 1/4, (1/512)^(1/4), 1/(2 * sqrt(8)), (1/2048)^(1/4), 1/8, (1/8192)^(1/4), 1/sqrt(128), 0];
P = quantize(W, codebook);
p1 = cat(1, repmat(P(1),numberofBeams,1), repmat(P(2),numberofBeams,1));
p_WB = diag(p1);
end

function ret = quantize(samples, codes)
    codes = codes(:);
    codes = sort(codes)';
    midpoints = (codes(2 : end) + codes(1 : end-1)) / 2;
    diffs = codes(2 : end) - codes(1 : end-1);
    ret = sum((samples(:) > midpoints) .* diffs, 2) + codes(1);
    ret = reshape(ret, size(samples));
end
