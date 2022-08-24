clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1层子带幅度加强的初始报告
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N1 = 8;                                                                    %水平方向天线数
N2 = 4;                                                                    %垂直方向天线数
v = 1;                                                                     %RI的值，由高层参数typeII-RI-Restriction-r16配置，UE不能报告RI > 4
Nt = N1 * N2;                                                              %天线数
P_CSIRS = 2 * N1 * N2;                                                     %端口数
phaseAlphabetSize = 4;                                                     %子带相位量化
subbandAmplitude = true;
N_SB = 24;                                                                 %传输带宽/子载波间隔
usernum = 10;                                                              %用户数
R = 1;                                                                     %高层参数numberofPMI-SubbandsPerCQI_Subband
paramCombination_r16 = 3;                                                  %高层R16组合参数，用于控制L,p_v和\beta
N3 = N_SB * R;
O3 = 4;
%选择组合参数和频域正交基
[L, p_v, beta, M_v, W_FD] = Combinationparam_and_FrequencyDomainBasis(paramCombination_r16, v, R, N3, O3);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sigma2 = 1;                                                                %信号功率????????????????????????????????????
Nr = 1;                                                                    %接收天线
Ncl = 5;                                                                   %簇数
Nray = 10;                                                                 %每簇的射线数                                                                
c = 3e8;
fc = 18e8;
lambda = c / fc;                                                           %信号波长
d = lambda / 2;                                                            %天线间距
M = 3;                                                                     %QAM调制阶数
snrdb = -5 : 1 : 20;
BER_vector = zeros(1, 26);
BER_vector1 = zeros(1, 26);
MSE = zeros(1,26);
N_3 = N_SB * R;


Length = 100;
loopnum = 1000;


% parfor (index = 1 : length(BER_vector), 12)
for index = 1 : length(BER_vector)
    for i = 1 : loopnum
        for loop1 = 1 : N_SB
            for loop2 = 1 : usernum
                  SNR = 10^(snrdb(index) / 10);
                  S = randi([0, 2^M - 1], N_SB, Length);
                  X = qammod(S, 2^M);
%生成S-V宽带信道
                  [H, H_WB, ~, ~] = wideband_mmwave_channel(N_SB, Nt, Nr, N1, N2, Ncl, Nray, lambda, sigma2, d, pi/4, pi/2);    %pi/4与pi/2为极化相位
                  power_H = H_WB * H_WB';
%选择波束（即空域压缩）并进行宽带幅度量化
                  [B1, B1_polar, ~, ~, W2_WB, b_max, b_max_polar] = W1_form_refine(L, N1, N2, Nr, Nt, H_WB);                    %获得最佳的2L个波束组
                  W1 = blkdiag(B1, B1_polar);        
                  W = zeros(2 * Nt, 1, N_SB);
                  Y = zeros(1, Length, N_SB);
                  Y1 = zeros(1, Length, N_SB);
                  W2_SB_FD = zeros(2 * L, N3);
                  d = zeros(1, O3);
                  d1 = zeros(1, N_SB);
%获取空频系数矩阵\widetilde{{W}_{2}^{SB}}
                  for k = 1 : N_SB
                      H_SB = H(:, :, k);
                      power_HSB = H_SB * H_SB';
%                       H_SB = H_SB / norm(H_SB, 'fro');
%                       W2_SB_FD(:, k) = pinv(W1 * W2_WB) * H_SB';
                      W2_SB_FD(:, k) = pinv(W1) * H_SB';
                  end
%从O3组正交基中选取一组
                  for i3 = 1 : O3
                     tempW = W2_SB_FD * W_FD(:, :, i3);
                     [~, D, ~] = svd(tempW);
                     d(i3) = D(1);
                  end
                  [~, O3pick] = max(d);
                  W_FD_i3 = W_FD(:, :, O3pick);
%选择W_f与生成\hat{{W}_{2}^{SB}}
                  for k = 1 : N_SB
                      temp_W = W2_SB_FD * W_FD_i3(:, k);
                      [~, D1, ~] = svd(temp_W);
                      d1(k) = D1(1);
                  end
                  [~, W_fpick] = maxk(d1, M_v);
                  W_f = W_FD_i3(:, W_fpick);
                  W2_SB_FDcompression = W2_SB_FD * W_f;
                  
%子带量化与码字生成
                  for k = 1 : N_SB
                      H_SB = H(:, :, k);      
                      p_SB = narrowband_quantization(W2_SB_FDcompression((1 : L), :), b_max);
                      p_SB_polar = narrowband_quantization(W2_SB_FDcompression((L + 1 : 2 * L), :), b_max_polar);
                      
%                       W2_SB = pinv(B1 * p_WB) * H_SB(:, (1 : Nt))';
%                       W2_SB_polar = pinv(B1_polar * p_WB_polar) * H_SB(:, (Nt + 1 : 2 * Nt))';
%                       W2_SB_FDbasis = blkdiag(W2_SB, W2_SB_polar);
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
                      Wmmse = (H1' * H1 + eye(1) / SNR) \ H1';
                      Wmmse = Wmmse / norm(Wmmse, 'fro');
                      Wmmse1 = (H11' * H11 + eye(1) / SNR) \ H11';
                      Wmmse1 = Wmmse1 / norm(Wmmse1, 'fro');
                      x = X(:, :, k);
%                       fan_x = norm(x,'fro');
%                       n = (1 / sqrt(2 * snr)) * (randn(1, Length) + 1i * randn(1, Length));
                      mse = sum((Wsb' - H_SB) .* conj(Wsb' - H_SB),'all') / numel(H_SB);
%                       y = Wmmse * H_SB * Wsb * x + n;  
                      y = awgn(Wmmse * H1 * x, SNR);
%                       y1 = Wmmse * H_SB * H_SB' * x + n;
                      y1 = awgn(Wmmse1 * H11 * x, SNR);
%                       y = H_SB * Wsb * x + n;                                
%                       y1 = H_SB * Wsb * x + n;
%                       fan_mmse=norm(Wmmse,'fro');
%                       fan_sb=norm(H_SB,'fro');
%                       fan_W=norm(Wsb,'fro');
                      Y(:, :, k) = y;
                      Y1(:,:,k) = y1;
                  end
            end
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
BER_vector = BER_vector / (Length * loopnum * N_SB);
BER_vector1 = BER_vector1 / (Length * loopnum * N_SB);
MSE = MSE / (loopnum * N_SB);
figure(1)
semilogy(snrdb, BER_vector, 'b-o',snrdb, BER_vector1, 'g-d');
legend('EtypeII码本','MF');
xlabel('SNR(dB)');
ylabel('BER');
hold on

