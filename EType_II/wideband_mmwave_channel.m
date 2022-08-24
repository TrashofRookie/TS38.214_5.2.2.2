function [H, H_WB, Ar, At] = wideband_mmwave_channel(L, Nt, Nr, N1, N2, Ncl, Nray, lambda, sigma2, d, gamma, eta)
%lambda                      波长
%theta                       方位角
%phi                         仰角
%sigma                       AOA/AOD的角度扩展
%d                           天线间距
%sigma2                      每簇的平均能量
% Initialization
mean_phi_aod = 0;
mean_theta_aoa = 0;
sigma_phi_aod = pi*50/180;
sigma_phi_aoa = pi*50/180;
sigma_theta_aod = pi*50/180;
H = zeros(Nr, 2 * Nt, L);
H_WB = zeros(Nr, 2 * Nt);

Ar = zeros(Nr, Ncl * Nray, L);
At = zeros(2 * Nt, Ncl * Nray, L);
phi_t = zeros(2 * Ncl * Nray, L);
theta_t = zeros(2 * Ncl * Nray, L);
theta_r = zeros(2 * Ncl * Nray, L);
P = [cos(gamma), exp(1i * eta) * sin(gamma)]';

  for l = 1 : L
      Hl = zeros(Nr, 2 * Nt);
%       Hl_with_beamSquint = zeros(Nr, Nt);
      index = 1;
      for tap = 1 : (Ncl * Nray)

          rayleigh_coeff = sqrt(sigma2/2)*(randn(1)+1j*randn(1));
          phi_t(tap, l) = genLaplacianSamples1(sigma_phi_aod);
          theta_r(tap, l) = genLaplacianSamples2(sigma_phi_aoa);
          theta_t(tap, l) = genLaplacianSamples2(sigma_theta_aod); 
          Ar(:, index, l) = angle(theta_r(tap, l), mean_theta_aoa, Nr, lambda, d);
          [at, ~, ~, ~, ~] = steervector(phi_t(tap, l), mean_phi_aod, theta_t(tap, l), lambda, Nt, N1, N2, d); 
          At(:, index, l) = kron(at, eye(2)) * P;
          Hl = Hl + rayleigh_coeff * Ar(:, index, l) * At(:, index, l)';
%               phi_t(tap, ray, l) = genLaplacianSamples1(sigma_phi_aod);
%               theta_r(tap, ray, l) = genLaplacianSamples2(sigma_phi_aoa);
%               theta_t(tap, ray, l) = genLaplacianSamples2(sigma_theta_aod);           
%               Ar(:, index, l) = kron(angle(theta_r(tap, ray, l), mean_theta_aoa, Nr, lambda, d), eye(2)) * P;
%               [at, ~, ~, ~, ~] = steervector(phi_t(tap, ray, l), mean_phi_aod, theta_t(tap, ray, l), lambda, Nt, N1, N2, d); 
%               At(:, index, l) = kron(at, eye(2)) * P;  
%               Hl = Hl + rayleigh_coeff * Ar(:, index, l) * At(:, index, l)';
          index = index + 1;
          H(:,:,l) = H(:,:,l) + Hl;

      end

      H(:,:,l) = sqrt((Nt * Nr) / (Nray * Ncl)) * H(:, :, l);     
      H_WB = H_WB + H(:, :, l);

  end

end

%生成收端导向矢量
function vectors_of_angles=angle(phi, phi0, M, lambda, array_element_spacing)

    %生成ULA相移                     
    wavenumber = 2 * pi / lambda; % k = 2pi/lambda
    phase_shift = wavenumber * array_element_spacing * sin(phi0 - phi) * (0 : M-1).';
    vectors_of_angles = exp(-1i * phase_shift);
end


% 基于逆变换采样生成随机变量
function x = genLaplacianSamples2(sigma_phi)
    u = rand(1,1);
    % 截断拉普拉斯的逆变换, \phi \in [-pi,pi]
    %sigma_phi   power azimuth spectrum (PAS)的标准偏移相位
    beta = 1/(1-exp(-sqrt(2)*pi/sigma_phi));
    x = beta*(exp(-sqrt(2)/sigma_phi*pi) - cosh(u));
end
function x = genLaplacianSamples1(sigma_phi)
    u = rand(1,1);
    % 截断拉普拉斯的逆变换, \phi \in [-pi,pi]
    %sigma_phi   power azimuth spectrum (PAS)的标准偏移相位
    beta = 1/(1-exp(-pi / (sqrt(2) * sigma_phi)));
    x = beta * (exp(-sqrt(2) / sigma_phi * pi) - cosh(u));
end


