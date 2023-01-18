function [mu_z, sigma_z] = log_normal_I_approx(Pt, Pi, dL, di, K, sigma_N, noise_dB, n)
%Fenton-Wilkinson Approximation for log-normal RV SINR
%Refer to I. Hadj-Kacem, H. Braham and S. B. Jemaa, "SINR and Rate 
%Distributions for Downlink Cellular Networks," in IEEE Transactions on 
%Wireless Communications, vol. 19, no. 7, pp. 4604-4616, July 2020.
%Date: 27/05/2021
%Author: Reuben Lim
%Date Modified: 9/6/2021
%Modified for: To account for Rician fading power as non central chi square
%To be paired with Gauss-Hermite quadrature method
%Approximates SINR as lognormal without the numerator Rician/Rayleigh
%fading. It is assumed that interference experiences Rician fading

%Pt = Tx power of SOI, Pi = Tx power of interferers
%dL = distance from GCS, di = vector of distances btw UAVs
%K = Rician factor, should be >= 5, n = path loss exponent (default 2)
%noise = noise power in dB (default -86)
if nargin < 8
    n = 2;
elseif nargin < 7
    n = 2;
    noise_dB = -86;
elseif nargin < 6
    n = 2;
    noise_dB = -86;
    sigma_N = 1;
end

hL = friis_cnst(Pt);
hi = friis_cnst(Pi);
noise = 10^(noise_dB/10)/1000;
eta = log(10)/10;
m = (gamma(1+1)/(1+K))*hypergeom(-1,1,-K);
var = (gamma(1+2)/(1+K)^2)*hypergeom(-2,1,-K) - m^2;
sigma = sqrt(var);
mu_ln = -(1/(2*eta))*log(sigma^2/m^4+1/m^2);
s_i = sqrt((2/eta^2)*(log(m)-eta*mu_ln)+sigma_N^2);
s_L = sigma_N;
v = mu_ln;
%Fenton-Wilkinson Method
theta_i = hi/hL; theta_0 = noise/hL;
temp_mz_i0 = exp(eta^2*s_L^2/2)*theta_0*dL^n;
temp_mz_in0 = exp(eta^2*(s_i^2+s_L^2)/2)*theta_i*sum((dL./di).^n,'all')*exp(eta*v);
m_z = temp_mz_i0 + temp_mz_in0;
%Finding the variance
temp1_i0 = exp(eta^2*s_L^2)*(exp(eta^2*s_L^2)-1);
temp1_in0 = exp(2*eta*v)*exp(eta^2*(s_i^2+s_L^2))*(exp(eta^2*(s_i^2+s_L^2))-1);
temp2_i0 = exp(1.5*eta^2*s_L^2)*(exp(eta^2*s_L^2)-1); %When i = 0
temp2_in0 = exp(eta^2*((s_i^2+s_L^2)/2+s_L^2))*(exp(eta^2*s_L^2)-1)*exp(eta*v); %When i not 0
temp3 = 0;
for i = 0:(length(di)-1)
    for j = (i+1):length(di)
        if i == 0
            temp3 = temp3 + theta_0*theta_i*(dL^2/di(j))^n*temp2_i0;
        else
            temp3 = temp3 + theta_i^2*(dL^2/(di(i)*di(j)))^n*temp2_in0;
        end
    end
end
temp3 = 2*temp3;
temp4 = temp1_i0*theta_0^2*dL^(2*n) + temp1_in0*theta_i^2*sum((dL./di).^(2*n),'all');
v_z = sqrt(temp3+temp4);
sigma_z = sqrt(log(v_z^2/m_z^2 + 1)/eta^2);
mu_z = (log(m_z)-(eta^2*sigma_z^2)/2)/eta;

end

function h_pl = friis_cnst(Pt)
    G_Tx = 1;
    G_Rx = 1;
    freq = 2.4e9;
    lambda = 3e8/freq;
    h_pl = Pt*(G_Tx*G_Rx*lambda^2)/(16*pi^2);
end