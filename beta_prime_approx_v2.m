function [k1,k2,theta1,theta2] = beta_prime_approx_v2(hL, hi, dL, di, K, omega, sigma_N, noise_dB, n)
%Beta prime approximation of SINR from ratio of Gamma distribution
%Date: 06/10/2021
%Author: Reuben Lim

%hL = Deterministic large-scale path loss of leader UAV, 
%hi = Deterministic large-scale path loss of member UAV, 
%dL = distance from GCS (single value)
%di = vector of distances btw UAVs (vector)
%sigma_N = logarithmic std dev of lognormal shadowing (vector, where the first element is sigma_N of SOI)
%K = Rician K factor of interference power (vector, where the first element is K of SOI)
%omega = Rician scale parameter (assuming same for all links) (single constant)
%n = path loss exponent (vector, where the first element is n of SOI)
%noise = noise power in dB (default -86)

%Checks that channel params vectors are correct:
num_I = length(di); %Number of interferers
assert((length(K) == num_I+1), "Rician K vector mismatch"); %Check Rician K factor vector
assert((length(sigma_N) == num_I+1), "sigma_N mismatch"); %Check sigma_N vector
assert((length(n) == num_I+1), "Path loss exponent vector mismatch"); %Check n vector

noise = 10^(noise_dB/10)/1000;
eta = log(10)/10;
mu_N = 0; % Logarithmic mean of lognormal shadowing

%Obtaining the mean and variance of interference
E_I = 0; % Mean of interference power received
var_I = 0; % Variance of interference power received
for i = 1:num_I %For each distance in interferer list
    mu_psi = mu_N/eta + 10*log10(hi*di(i)^(-n(i+1))); %Logarithmic mean of large-scale power from ith interferer
    sigma_psi = sigma_N(i+1); %Logarithmic std dev of large-scale power from ith interferer
    E_psi = exp(eta*mu_psi + eta^2*sigma_psi^2/2); %Mean of large-scale power from ith interferer
    var_psi = exp(2*eta*mu_psi+eta^2*sigma_psi^2)*(exp(eta^2*sigma_psi^2)-1); %Variance of large-scale power from ith interferer
    Ki = K(i+1);
    E_chi = (gamma(1+1)/(1+Ki))*hypergeom(-1,1,-Ki)*omega; %Mean of Rician fading for ith interferer
    var_chi = (gamma(1+2)/(1+Ki)^2)*hypergeom(-2,1,-Ki)*omega^2 - E_chi^2; %Variance of Rician fading for ith interferer
    E_I = E_I + E_psi*E_chi;
    var_I = var_I + (var_psi+E_psi^2)*(var_chi+E_chi^2) - E_psi^2*E_chi^2;
end
%Approximating the denominator of SINR as a gamma dist.
E_X_den = E_I + noise; % Mean of denominator of SINR
var_X_den = var_I; % Variance of denominator of SINR
k2 = E_X_den^2/var_X_den; % Shape factor for Gamma approx. of denominator of SINR
theta2 = var_X_den/E_X_den; % Scale factor for Gamma approx. of denominator of SINR

%Approximating the numerator of SINR as a gamma dist.
mu_psi_soi = mu_N + 10*log10(hL*dL^(-n(1))); %Logarithmic mean of large-scale power from SOI
sigma_psi_soi = sigma_N(1); %Logarithmic std dev of large-scale power from SOI
E_chi_soi = (gamma(1+1)/(1+K(1)))*hypergeom(-1,1,-K(1))*omega; %Mean of Rician fading for SOI
var_chi_soi = (gamma(1+2)/(1+K(1))^2)*hypergeom(-2,1,-K(1))*omega^2 - E_chi_soi^2; %Variance of Rician fading for SOI
E_psi_soi = exp(eta*mu_psi_soi + eta^2*sigma_psi_soi^2/2); % Mean of shadowing in numerator
var_psi_soi = exp(2*eta*mu_psi_soi+eta^2*sigma_psi_soi^2)*(exp(eta^2*sigma_psi_soi^2)-1); % Variance of shadowng in numerator
E_X_num = E_psi_soi*E_chi_soi; % Mean of numerator of SOI
var_X_num = (var_psi_soi+E_psi_soi^2)*(var_chi_soi+E_chi_soi^2) - E_psi_soi^2*E_chi_soi^2; % Variance of numerator of SOI
k1 = E_X_num^2/var_X_num; % Shape factor for Gamma approx. of denominator of SINR
theta1 = var_X_num/E_X_num; % Scale factor for Gamma approx. of denominator of SINR

end

% function h_pl = friis_cnst(Pt)
%     G_Tx = 1;
%     G_Rx = 1;
%     freq = 2.4e9;
%     lambda = 3e8/freq;
%     h_pl = Pt*(G_Tx*G_Rx*lambda^2)/(16*pi^2);
% end