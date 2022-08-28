clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1层子带幅度加强的初始报告
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N1 = 8;
N2 = 4;
RI = 1;                                                                    %UE不能报告RI>2
Nt = N1*N2;
P_CSIRS = 2 * N1 * N2;
numberofBeams = 4;
phaseAlphabetSize = 4;
subbandAmplitude = true;
numberofSubbands = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sigma2 = 1;                                                                %信号功率
Nr = 1;                                                                    %接收天线
Ncl = 5;                                                                   %簇数
Nray = 5;                                                                  %每簇的射线数                                                                
c = 3e8;
fc = 18e8;
lambda = c / fc;                                                           %信号波长
d = lambda / 2;                                                            %天线间距
M = 6;                                                                     %QAM调制阶数
snrdb = -5 : 1 : 20;
BER_vector = zeros(1, 26);
BER_vector1 = zeros(1, 26);
BER_vector11 = zeros(1, 26);
MSE = zeros(1,26);
Length = 1000;
loopnum = 1000;
parfor (index = 1 : length(BER_vector), 12)
%for index = 1 : length(BER_vector)
    for i = 1 : loopnum
        snr = 10^(snrdb(index) / 10);
        S = randi([0, 2^M - 1], 1, Length, numberofSubbands);
        X = qammod(S, 2^M);
        [H, H_WB, ~, ~] = wideband_mmwave_channel(numberofSubbands, Nt, Nr, N1, N2, Ncl, Nray, lambda, sigma2, d, pi/4, pi/2);     %pi/4与pi/2为极化相位
        [B, W1, p_WB, beam_max, ~] = beamselect_and_W1form(numberofBeams, N1, N2, Nr, Nt, H_WB);                                   %获得最佳的2L个波束组    
        W = zeros(2 * Nt, 1, numberofSubbands);
        Y = zeros(1, Length, numberofSubbands);
        Y1 = zeros(1, Length, numberofSubbands);
        Y11 = zeros(1, Length, numberofSubbands);
        for k = 1 : numberofSubbands
            H_SB = H(:, :, k);
            Rn = H_SB' * H_SB;
            hat_Rn = W1' * Rn * W1;
            [optimal_BCC, ~] = eig(hat_Rn);
            optimal_BCC = optimal_BCC(:,1);
            optimal_BCC = optimal_BCC / norm(optimal_BCC, 'fro');
            H_SB = H_SB / norm(H_SB, 'fro');
            [~, ~, V] = svd(H_SB);
            V_SB = V(:, 1);                                                %子带反馈矩阵
            W2_SB = W1' * V_SB;
            p_SB = narrowband_quantization(W2_SB);
            c = narrowband_phasequan(W2_SB, phaseAlphabetSize);
            gamma = N1 * N2 * sum((p_WB * p_SB).^2, 'all');
            gamma = 1 / sqrt(gamma);
            W(:, :, k) = gamma * W1 * p_WB * p_SB * c;
            W_optimal_BCC = W1 * optimal_BCC;
            Wsb = W(:, :, k);
            Wsb = Wsb/norm(Wsb, 'fro');
            H1 = H_SB * Wsb;
            H1 = H1/norm(H1,'fro');
            H11 = H_SB * H_SB';
            H11 = H11/norm(H11,'fro');
            H111 = H_SB * W_optimal_BCC;
            H111 = H111/norm(H111,'fro');
            Wmmse = (H1' * H1 + eye(1) / snr) \ H1';
            Wmmse = Wmmse / norm(Wmmse, 'fro');
            Wmmse1 = (H11' * H11 + eye(1) / snr) \ H11';
            Wmmse1 = Wmmse1 / norm(Wmmse1, 'fro');
            Wmmse11 = (H111' * H111 + eye(1) / snr) \ H111';
            Wmmse11 = Wmmse11 / norm(Wmmse11, 'fro');
            x = X(:, :, k);
            fan_x = norm(x,'fro');
            mse = sum((Wsb' - V_SB) .* conj(Wsb' - V_SB),'all') / numel(V_SB);
            y = awgn(Wmmse * H1 * x, snr);
            y1 = awgn(Wmmse1 * H11 * x, snr);
            %y11 = awgn(Wmmse11 * H111 * x, snr);
            fan_mmse = norm(Wmmse,'fro');
            fan_sb = norm(H_SB,'fro');
            fan_W =norm(Wsb,'fro');
            fan_VSB = norm(V_SB,'fro');
            Y(:, :, k) = y;
            Y1(:, :, k) = y1;
            %Y11(:, :, k) = y11;
        end
        receive = qamdemod(Y, 2^M);
        receive1 = qamdemod(Y1,2^M);
        %receive11 = qamdemod(Y11,2^M);
        Sum = sum(receive ~= S, 'all');
        Sum1= sum(receive1 ~= S, 'all');
        %Sum11= sum(receive11 ~= S, 'all');
        BER_vector(index) = BER_vector(index) + Sum;
        BER_vector1(index) = BER_vector1(index) + Sum1;
        %BER_vector11(index) = BER_vector11(index) + Sum11;
    end
end
BER_vector = BER_vector / (Length * loopnum * numberofSubbands);
BER_vector1 = BER_vector1 / (Length * loopnum * numberofSubbands);
%BER_vector11 = BER_vector11 / (Length * loopnum * numberofSubbands);
figure(1)
semilogy(snrdb, BER_vector, 'b-o',snrdb, BER_vector1, 'g-d');
legend('TypeII码本','MF','BCC');
xlabel('SNR(dB)');
ylabel('BER');
hold on
