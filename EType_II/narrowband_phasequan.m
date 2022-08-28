function [C] = narrowband_phasequan(W, N_PSK)
W_phase = angle(W);
% W_phase_angle = W_phase*180/pi;
% W_phase = exp(1i .* W_phase);
if N_PSK == 4
    codebook = pi * (-1 : 2) /2;
elseif N_PSK == 8
    codebook = pi * (-3 : 4) / 4;
end
C = exp(1i .* quantize(W_phase, codebook));
end

function ret = quantize(samples, codes)
    codes = codes(:);
    codes = sort(codes)';
    midpoints = (codes(2 : end) + codes(1 : end-1)) / 2;
    diffs = codes(2 : end) - codes(1 : end-1);
    ret = sum((samples(:) > midpoints) .* diffs, 2) + codes(1);
    ret = reshape(ret, size(samples));
end