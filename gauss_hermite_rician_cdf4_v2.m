function sinr_cdf = gauss_hermite_rician_cdf4_v2(c, mu_ln, sigma_ln, K, omega, Np)
%Date: 28/09/2021
%Author: Reuben Lim
%Desc: Applying Gauss-Hermite Quadrature to approx. PDF of product 
%      distribution of Rician squared RV and lognormal RV for the points c. 
% c        : Points to evaluate PDF at (input points)
% mu_ln    : Logarithmic mean of lognormal RV
% sigma_ln : Logarithmic std dev of lognormal RV
% K        : Rician K factor of Rician squared RV
% omega    : Scale parameter of Rician squared RV
% Np       : Number of points to evaluate Gauss-Hermite Quadrature

if nargin < 6
    Np = 95;
elseif nargin < 5
    Np = 95;
    omega = 1;
end

%Gauss Hermite Approximation
[xi,wi] = GaussHermite_2(32); % Fix the Gauss Hermite terms at 32
sinr_cdf = zeros(1,length(c));
inf_trig = 0; %This flag will be set if Np has been found to be insufficient
for i = 1:length(c)
    %i
%     This might not be a good way to check for -inf cases
%     if Np <= 95
%         sinr_cdf(i) = F1(c(i),xi,wi,K,omega,mu_ln,sigma_ln,Np);
%     else
%         sinr_cdf(i) = py.gh_rician_cdf_F1.F1(c(i),xi,wi,K,omega,mu_ln,sigma_ln,Np);
%     end
    
    if inf_trig == 0
        sinr_cdf(i) = F1(c(i),xi,wi,K,omega,mu_ln,sigma_ln,Np);
    elseif inf_trig == 1
        %sprintf("Progress: %d / %d", i, length(c));
        sinr_cdf(i) = py.gh_rician_cdf_F1.F1(c(i),xi,wi,K,omega,mu_ln,sigma_ln,Np);
        %sinr_cdf(i) = F1_vpa(c(i),xi,wi,K,omega,mu_ln,sigma_ln,Np);
    end
    
    if isnan(sinr_cdf(i)) || isinf(sinr_cdf(i))
        %sprintf("Progress: %d / %d", i, length(c));
        sinr_cdf(i) = py.gh_rician_cdf_F1.F1(c(i),xi,wi,K,omega,mu_ln,sigma_ln,Np);
        %sinr_cdf(i) = F1_vpa(c(i),xi,wi,K,omega,mu_ln,sigma_ln,Np);
        inf_trig = 1;
    end
    %sinr_cdf(i) = py.gh_rician_cdf_F1.F1(c(i),xi,wi,K,omega,mu_ln,sigma_ln,Np);
    %sinr_cdf(i) = F1(c(i),xi,wi,K,omega,mu_ln,sigma_ln,Np);
end

end

function f = F1(x,xi,wi,K,omega,m_z,v_z,N)
    %Gauss-Hermite f(x) function for hybrid approx PDF with lognormal
    %x is a single values
    eta = log(10)/10;
    b = (K+1)/omega;
    
%     z = vpa(exp(sqrt(2)*eta*v_z.*xi + eta*m_z));
%     x = vpa(x);
    z = exp(sqrt(2)*eta*v_z.*xi + eta*m_z);
    temp = 0;
    for i = 0:N
        for j = 0:i
            E_ln = sum(wi./((z).^(j)).*exp(-b.*x./z))/sqrt(pi);
            %temp = temp + ((vpa(K)^i*vpa(b)^j)/(factorial(vpa(i))*factorial(vpa(j)))) * exp(-K) * x^j * E_ln;
            temp = temp + (K^i*b^j)/(factorial(i)*factorial(j)) * exp(-K) * x^j * E_ln;
        end
    end
    f = 1 - temp;
end

function f = F1_vpa(x,xi,wi,K,omega,m_z,v_z,N)
    % This is the vpa version
    %Gauss-Hermite f(x) function for hybrid approx PDF with lognormal
    %x is a single values
    eta = log(10)/10;
    b = (K+1)/omega;
    
    z = vpa(exp(sqrt(2)*eta*v_z.*xi + eta*m_z));
    x = vpa(x);
%     z = exp(sqrt(2)*eta*v_z.*xi + eta*m_z);
    temp = 0;
    for i = 0:N
        for j = 0:i
            E_ln = sum(wi./((z).^(j)).*exp(-b.*x./z))/sqrt(pi);
            temp = temp + ((vpa(K)^i*vpa(b)^j)/(factorial(vpa(i))*factorial(vpa(j)))) * exp(-K) * x^j * E_ln;
            %temp = temp + (K^i*b^j)/(factorial(i)*factorial(j)) * exp(-K) * x^j * E_ln;
        end
    end
    f = 1 - temp;
end