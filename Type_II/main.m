clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1层子带幅度加强的初始报告
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N1 = 8;
N2 = 2;
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
M = 4;                                                                     %QAM调制阶数
snrdb = -5 : 1 : 20;
BER_vector = zeros(1, 26);
BER_vector1 = zeros(1, 26);
MSE = zeros(1,26);
Length = 100;
loopnum = 1000;
parfor (index = 1 : length(BER_vector), 12)
%for index = 1 : length(BER_vector)
    for i = 1 : loopnum
        snr = 10^(snrdb(index) / 10);
        S = randi([0, 2^M - 1], 1, Length, numberofSubbands);
        X = qammod(S, 2^M);
        [H, H_WB, ~, ~] = wideband_mmwave_channel(numberofSubbands, Nt, Nr, N1, N2, Ncl, Nray, lambda, sigma2, d, pi/4, pi/2);                 %pi/4与pi/2为极化相位
        [B1, B1_polar, p_WB, p_WB_polar, p_WB_quan, p_WB_polar_quan, ~, ~] = W1_form_refine(numberofBeams, N1, N2, Nr, Nt, H_WB);%获得最佳的2L个波束组
        W1 = [B1, zeros(size(B1)); zeros(size(B1)), B1_polar];
        P_WB = [p_WB_quan, zeros(size(p_WB_quan)); zeros(size(p_WB_quan)), p_WB_polar_quan];         
        W = zeros(2 * Nt, 1, numberofSubbands);
        Y = zeros(1, Length, numberofSubbands);
        Y1 = zeros(1, Length, numberofSubbands);
        for k = 1 : numberofSubbands
            H_SB = H(:, :, k);      
            H_SB = H_SB / norm(H_SB, 'fro');
            W2_SB = pinv(B1) * H_SB(:, (1 : Nt))';
            W2_SB_polar = pinv(B1_polar) * H_SB(:, (Nt + 1 : 2 * Nt))';
            p_SB = narrowband_quantization(W2_SB);
            p_SB_polar = narrowband_quantization(W2_SB_polar);
            P_SB = [p_SB, zeros(size(p_SB)); zeros(size(p_SB)), p_SB_polar];
            c = narrowband_phasequan(W2_SB, phaseAlphabetSize);
            c_polar = narrowband_phasequan(W2_SB_polar, phaseAlphabetSize);
            C = cat(1, c, c_polar);
            gamma = N1 * N2 * sum((P_WB * P_SB).^2, 'all');
            gamma = 1 / sqrt(gamma);
            W(:, :, k) = gamma * W1 * P_WB * P_SB * C;
            Wsb = W(:, :, k);
            Wsb = Wsb/norm(Wsb, 'fro');
            H1 = H_SB * Wsb;
            H11 = H_SB * H_SB';
            Wmmse = (H1' * H1 + eye(1) / snr) \ H1';
            Wmmse = Wmmse / norm(Wmmse, 'fro');
            Wmmse1 = (H11' * H11 + eye(1) / snr) \ H11';
            Wmmse1 = Wmmse1 / norm(Wmmse1, 'fro');
            x = X(:, :, k);
%             fan_x = norm(x,'fro');
%             n = (1 / sqrt(2 * snr)) * (randn(1, Length) + 1i * randn(1, Length));
            mse = sum((Wsb' - H_SB) .* conj(Wsb' - H_SB),'all') / numel(H_SB);
%             y = Wmmse * H_SB * Wsb * x + n;  
            y = awgn(Wmmse * H1 * x, snr);
%             y1 = Wmmse * H_SB * H_SB' * x + n;
            y1 = awgn(Wmmse1 * H11 * x, snr);
%             y = H_SB * Wsb * x + n;                                
%             y1 = H_SB * Wsb * x + n;
%             fan_mmse=norm(Wmmse,'fro');
%             fan_sb=norm(H_SB,'fro');
%             fan_W=norm(Wsb,'fro');
            Y(:, :, k) = y;
            Y1(:,:,k) = y1;
        end
        receive = qamdemod(Y, 2^M);
        receive1 = qamdemod(Y1,2^M);
        Sum = sum(receive ~= S, 'all');
        Sum1= sum(receive1 ~= S, 'all');
        BER_vector(index) = BER_vector(index) + Sum;
        BER_vector1(index) = BER_vector1(index) + Sum1;
        MSE(index) = MSE(index) + mse;
    end
end
BER_vector = BER_vector / (Length * loopnum * numberofSubbands);
BER_vector1 = BER_vector1 / (Length * loopnum * numberofSubbands);
MSE = MSE / (loopnum * numberofSubbands);
figure(1)
semilogy(snrdb, BER_vector, 'b-o',snrdb, BER_vector1, 'g-d');
legend('TypeII码本','MF');
xlabel('SNR(dB)');
ylabel('BER');
hold on
