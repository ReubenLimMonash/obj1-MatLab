function [mu_psi_j, sigma_psi_j] = log_normal_large_scale_approx(hL, hi, dL, di, K, omega, sigma_N, noise_dB, n)
%Date: 21/09/2021
%Desc: Moment matching approach to approximate interference and large-scale
%      power received as lognormal distribution.

%hL = Deterministic large-scale path loss of leader UAV, 
%hi = Deterministic large-scale path loss of member UAV, 
%dL = distance from GCS
%di = vector of distances btw UAVs
%sigma_N = logarithmic std dev of lognormal shadowing (vector, where the first element is sigma_N of SOI)
%K = Rician K factor of interference power (vector) 
%omega = Rician scale parameter (assuming same for all links) (single constant)
%n = path loss exponent (vector, where the first element is n of SOI)
%noise = noise power in dB (default -86)

%Checks that channel params vectors are correct:
num_I = length(di); %Number of interferers
assert((length(K) == num_I), "Rician K vector mismatch"); %Check Rician K factor vector
assert((length(sigma_N) == num_I+1), "sigma_N mismatch"); %Check sigma_N vector
assert((length(n) == num_I+1), "Path loss exponent vector mismatch"); %Check n vector

% hL = friis_cnst(Pt);
% hi = friis_cnst(Pi);
noise = 10^(noise_dB/10)/1000;
eta = log(10)/10;
mu_N = 0; % Logarithmic mean of lognormal shadowing

%Approximate the interference as lognormal dist.
E_I = 0; % Mean of interference power received
var_I = 0; % Variance of interference power received
for i = 1:num_I %For each distance in interferer list
    mu_psi = mu_N/eta + 10*log10(hi*di(i)^(-n(i+1))); %Logarithmic mean of large-scale power from ith interferer
    sigma_psi = sigma_N(i+1); %Logarithmic std dev of large-scale power from ith interferer
    E_psi = exp(eta*mu_psi + eta^2*sigma_psi^2/2); %Mean of large-scale power from ith interferer
    var_psi = exp(2*eta*mu_psi+eta^2*sigma_psi^2)*(exp(eta^2*sigma_psi^2)-1); %Variance of large-scale power from ith interferer
    Ki = K(i);
    E_chi = (gamma(1+1)/(1+Ki))*hypergeom(-1,1,-Ki)*omega; %Mean of Rician fading for ith interferer
    var_chi = (gamma(1+2)/(1+Ki)^2)*hypergeom(-2,1,-Ki)*omega^2 - E_chi^2; %Variance of Rician fading for ith interferer
    E_I = E_I + E_psi*E_chi;
    var_I = var_I + (var_psi+E_psi^2)*(var_chi+E_chi^2) - E_psi^2*E_chi^2;
end
sigma_I = sqrt(log(var_I/E_I^2 + 1)/eta^2);
mu_I = (log(E_I)-eta^2*sigma_I^2/2)/eta;

%Approximating the large-scale gain of received power and interference as
%lognormal
	% Parameters ----------------------------------------------------------
mu_psi_soi = mu_N + 10*log10(hL*dL^(-n(1))); %Logarithmic mean of large-scale power from SOI
sigma_psi_soi = sigma_N(1); %Logarithmic std dev of large-scale power from SOI
mu_psi_1 = 10*log10(noise) - mu_psi_soi; % Psi 1 accounting for noise
mu_psi_2 = mu_I - mu_psi_soi; %Psi 2 accounting for interference
sigma_psi_1 = sigma_psi_soi; 
sigma_psi_2 = sqrt(sigma_I^2 + sigma_psi_soi^2);
    % ---------------------------------------------------------------------
E_psi_1 = exp(eta*mu_psi_1 + eta^2*sigma_psi_1^2/2); % Mean of Psi 1
E_psi_2 = exp(eta*mu_psi_2 + eta^2*sigma_psi_2^2/2); % Mean of Psi 2
var_psi_1 = exp(2*eta*mu_psi_1+eta^2*sigma_psi_1^2)*(exp(eta^2*sigma_psi_1^2)-1); %Variance of Psi 1
var_psi_2 = exp(2*eta*mu_psi_2+eta^2*sigma_psi_2^2)*(exp(eta^2*sigma_psi_2^2)-1); %Variance of Psi 2
E_psi_j = E_psi_1 + E_psi_2; % Mean of "overall" Psi
cov_psi_1_2 = exp(eta*(mu_psi_1+mu_psi_2) + eta^2*(sigma_psi_1^2+sigma_psi_2^2+2*sigma_psi_soi^2)/2)...
                - E_psi_1 * E_psi_2; % Covariance pf Psi 1 and Psi 2
var_psi_j = var_psi_1 + var_psi_2 + 2*cov_psi_1_2; % Variance of "overall" Psi

sigma_psi_j = sqrt(log(var_psi_j/E_psi_j^2+1)/eta^2);% Logarithmic variance of "overall" Psi
mu_psi_j = -(log(E_psi_j)-eta^2*sigma_psi_j^2/2)/eta; % Logarithmic mean of "overall" Psi
%NOTE: mu_psi_j has already been multiplied by -1 in function.
end

% function h_pl = friis_cnst(Pt)
%     G_Tx = 1;
%     G_Rx = 1;
%     freq = 2.4e9;
%     lambda = 3e8/freq;
%     h_pl = Pt*(G_Tx*G_Rx*lambda^2)/(16*pi^2);
% end