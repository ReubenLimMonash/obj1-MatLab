function sinr_pdf = gauss_hermite_rician_pdf_v2(c, mu_ln, sigma_ln, K, omega, Np)
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
    Np = 32;
elseif nargin < 5
    Np = 32;
    omega = 1;
end

%Gauss Hermite Approximation
[xi,wi] = GaussHermite_2(Np);
sinr_pdf = zeros(1,length(c));
for i = 1:length(c)
    f = f1(c(i),xi,K,omega,mu_ln,sigma_ln);
    sinr_pdf(i) = sum(wi.*f)/sqrt(pi);
end

end

function y = rician_pwr_pdf(x,K,omega)
    x = vpa(x);
    y = (K+1)/omega .* exp(-K-(K+1).*x./omega).*besseli(0, 2*sqrt((K^2+K).*x./omega));
end

function y = rician_pwr_pdf_approx(x,K,omega)
    %Using the sum of scaled exponential approx. for besseli
    x = vpa(x);
    [a,b] = besseli0_approx(x);
    besseli_approx = 0;
    for i = 1:4
        besseli_approx = besseli_approx + a(i)*exp(b(i)*2*sqrt((K^2+K).*x./omega));
    end
    y = (K+1)/omega .* exp(-K-(K+1).*x./omega).*besseli_approx;
end

function f = f1(x,xi,K,omega,m_z,v_z)
    %Gauss-Hermite f(x) function for hybrid approx PDF with lognormal
    eta = log(10)/10;
    z = exp(sqrt(2)*eta*v_z.*xi + eta*m_z);
    y = rician_pwr_pdf(x./z,K,omega);
    %y = rician_pwr_pdf_approx(x./z,K1,1);
    f = (1./z) .* y;
end