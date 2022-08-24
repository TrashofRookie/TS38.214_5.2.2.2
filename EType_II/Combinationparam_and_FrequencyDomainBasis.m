function [L, p_v, beta, M_v, W_FD] = Combinationparam_and_FrequencyDomainBasis(paramCombination_r16, v, R, N3, O3)
table5_2_2_2_5_1 = [
    2 1/4 1/8 1/4;
    2 1/4 1/8 1/2;
    4 1/4 1/8 1/4;
    4 1/4 1/8 1/2;
    4 1/4 1/4 3/4;
    4 1/2 1/4 1/2;
    6 1/4 0 1/2;
    6 1/4 0 3/4];
L = table5_2_2_2_5_1(paramCombination_r16, 1);
p_v = (v<=2) .* table5_2_2_2_5_1(paramCombination_r16, 2) + (v>2) .* table5_2_2_2_5_1(paramCombination_r16, 3);
beta = table5_2_2_2_5_1(paramCombination_r16, 4);
M_v = ceil(p_v * N3 / R);
W_FD = zeros(N3, N3, O3);
n3 = 0 : 1 : (N3 - 1);
for q = 1 : O3
        W_FD(:, :, q) = exp(1i * 2 * pi *((n3' *O3 + (q - 1)) * n3) / (N3 * O3));
end
end