function [mu, sigma] = lognormal_rice_moment_match(mu_ln, sigma_ln, K, omega)
    % This function approximates the product distribution of a lognormal
    % distribution and a Rician squared distribution using a lognormal
    % distribution by moment matching. 
    % mu_ln    : Logarithmic mean of lognormal RV (single value)
    % sigma_ln : Logarithmic std dev of lognormal RV (single value)
    % K        : Rician K factor of Rician squared RV (single value)
    % omega    : Scale parameter of Rician squared RV (single value)
    %
    % Returns:-
    % mu       : Logarithmic mean of lognormal approximation
    % sigma    : Logarithmic std dev of lognormal approximation
    
    E_ln = exp(mu_ln + sigma_ln^2/2); %Mean of LN RV
    var_ln = exp(2*mu_ln+sigma_ln^2)*(exp(sigma_ln^2)-1); %Variance of LN RV
    E_chi = (gamma(1+1)/(1+K))*hypergeom(-1,1,-K)*omega; %Mean of Rician squared RV
    var_chi = (gamma(1+2)/(1+K)^2)*hypergeom(-2,1,-K)*omega^2 - E_chi^2; %Variance of Rician squared RV
    E_x = E_ln * E_chi;
    var_x = (var_ln+E_ln^2)*(var_chi+E_chi^2) - E_ln^2*E_chi^2;
     
    eta = log(10)/10;
    sigma = sqrt(log(var_x/E_x^2 + 1)/eta^2);
    mu = (log(E_x)-eta^2*sigma^2/2)/eta;
end