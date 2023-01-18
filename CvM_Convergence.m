%Date: 23/7/2021
%Testing for number of histogram bins needed for convergence of CvM measure
%of CVM vs K and sigma plot (Fig. 1).

%Clear and load the Python func
clear classes;
pyfunc = py.importlib.import_module('gh_rician_cdf_F1');
py.importlib.reload(pyfunc);

% Specify the input parameters ===========================================
Pt_dBm = 40; Pt = 10^(Pt_dBm/10)/1000; %40dBm = 10W Tx power
Pi_dBm = 10; Pi = 10^(Pi_dBm/10)/1000; %23dBm = 0.1995W Tx power
freq = 0.968e9; %Channel frequency
height = 100; %Height btw UAVs and GCS (m)
dL_H = 50; %Horizontal distance between leader UAV and GCS
dL = sqrt(dL_H^2 + height^2);
di = [5 5 5]; %Distances between each follower UAV and the leader UAV
noise_dB = -107;
noise = 10^(noise_dB/10)/1000;
omega = 1;   %Rician Scale factor
theta_GCS = 45; % Antenna tilt angle of GCS (in degrees)
Go = 1; % Antenna gain at horizontal
%Environment parameters-----------------------------
%ENVIRONMENT = SUBURBAN
n_max = 2.75; n_min = 2; %Range of path loss exponent
K_dB_max = 17.5; K_dB_min = 7.8; %Range of Rician K (dB) in suburban Cleveland
alpha = 11.25; beta = 0.06; %Env parameters for logarithm std dev of shadowing 
sigma_N_i_dB = -inf;
%Calculate channel parameters affected by positions
PLoS = PLoS_v3(dL_H,0,height,0.1,7.5e-4,8); % LoS Prob. btw leader & GCS
theta_L_d = atand(height/dL_H); %Elevation angle of leader from GCS in degrees
G_GCS_L = antenna_gain(theta_L_d,theta_GCS,Go); % Antenna gain from GCS to leader
G_L_GCS = antenna_gain(theta_L_d,0,Go); % Antenna gain from leader to GCS
G_i = antenna_gain(0,0,Go); % Antenna gain between leader and member i
hL = friis_cnst(Pt,G_GCS_L,G_L_GCS,freq); %Channel gain U2G
hi = friis_cnst(Pi,G_i,G_i,freq); %Channel gain U2U
n_L = (n_min-n_max)*PLoS + n_max; %Path loss exponent btw leader and GCS
n_i = (n_min-n_max)*1 + n_max; %Path loss exponent btw leader and member i (PLoS = 1)
%sigma_N_L_dB = alpha*exp(-beta*theta_L_d); sigma_N_L = 10^(sigma_N_L_dB/10);%Logarithmic Std dev of shadowing btw leader & GCS
sigma_N_i = 10^(sigma_N_i_dB/10);%Logarithmic Std dev of shadowing btw leader & member i
%K_L_dB = K_dB_min*exp(log(K_dB_max/K_dB_min)*PLoS^2); K_L = 10^(K_L_dB/10); %Rician K factor btw leader and GCS
K_i_dB = K_dB_max; K_i = 10^(K_i_dB/10); %Rician K factor btw leader and member i
eta = log(10)/10;
N = 10^8; %Number of Monte Carlo data points
%Np = 32; %Number of Gauss-Hermite points to use

%Custom shadowing / multipath fading
sigma_N_L_dB = 0; sigma_N_L = 10^(sigma_N_L_dB/10);
%K_L_dB = 14; K_L = 10^(K_L_dB/10); 

num_points = [50, 75, 100, 125, 150, 175];
K_L_dB = [8, 8, 16, 16];
sigma_N_L_dB = [0, 8, 0, 8];
CVM_ln = zeros(length(K_L_dB),length(sigma_N_L_dB),length(num_points)); % Cramer-von Mises criterion for LN
CVM_beta = zeros(length(K_L_dB),length(sigma_N_L_dB),length(num_points)); % Cramer-von Mises criterion for BP
CVM_gh = zeros(length(K_L_dB),length(sigma_N_L_dB),length(num_points)); % Cramer-von Mises criterion for GH
for k = 1:length(K_L_dB)
    for j = 1:length(sigma_N_L_dB)
        K_L = 10^(K_L_dB(k)/10); 
        sigma_N_L = 10^(sigma_N_L_dB(j)/10); 
        sprintf("K = %d dB, \x3C3 = %d dB",K_L_dB(k),sigma_N_L_dB(j));
        v_L = sqrt(omega*K_L/(1+K_L)); s_L = sqrt(omega/(2+2*K_L)); %Rician parameters for leader UAV
        v_i = sqrt(omega*K_i/(1+K_i)); s_i = sqrt(omega/(2+2*K_i)); %Rician parameters for member UAV
        Pwr_sim = hL*dL^(-n_L).*ricePwrrnd(v_L*ones(1, N),s_L).*10.^(normrnd(0,sigma_N_L,1,N)./10);
        I_sim = ones(length(di), N);
        for i = 1:length(di)
            I_sim(i,:) = (di(i)).^(-n_i).*ricePwrrnd(v_i*ones(1, N),s_i).*10.^(normrnd(0,sigma_N_i,1,N)./10);
        end
        I_sim = hi.*sum(I_sim);
        SINR_sim = Pwr_sim./(I_sim+noise); %SINR
        %Get simulated SINR CDF
        mean_sinr = mean(SINR_sim); %Mean of SINR
        sigma_sinr = std(SINR_sim); %Std dev of SINR
        % Finding convergence
        for i = 1:length(num_points)
            h = histogram(SINR_sim,num_points(i),'BinLimits',[max(0,mean_sinr-5*sigma_sinr),mean_sinr+5*sigma_sinr],'Normalization','cdf');
            CDF_sinr_sim = h.Values;
            c = conv(h.BinEdges, [0.5 0.5], 'valid'); %CDF points
            c_index = CDF_sinr_sim > 0;
            sinr_th = c(c_index); 
            CDF_sinr_sim = CDF_sinr_sim(c_index);
            %======================================================================
            % Approximations of SINR
            K = [K_L, K_i*ones(1,length(di))]; %Corresponding to [GCS, Member 1, Member 2]
            sigma_N = [sigma_N_L, sigma_N_i*ones(1,length(di))]; %Corresponding to [GCS, Member 1, Member 2]
            n = [n_L, n_i*ones(1,length(di))]; %Corresponding to [GCS, Member 1, Member 2]
            % Get the lognormal approximation CDF
            [mu_z_0,sigma_z_0] = log_normal_large_scale_approx(hL, hi, dL, di, K(2:end), omega, sigma_N, noise_dB, n);
            [mu_z_ln,sigma_z_ln] = lognormal_rice_moment_match(eta*mu_z_0,eta*sigma_z_0,K(1),omega);
            CDF_sinr_ln = log_normal_cdf(sinr_th, eta*mu_z_ln, eta*sigma_z_ln);
            CDF_sinr_ln(isnan(CDF_sinr_ln))=0;
            err_cdf_ln = (CDF_sinr_ln-CDF_sinr_sim);
            err_cdf_ln(isnan(err_cdf_ln))=CDF_sinr_ln(isnan(err_cdf_ln));
            err_cdf_ln(isinf(err_cdf_ln))=CDF_sinr_ln(isinf(err_cdf_ln));
            CVM_ln(k,j,i) = sum((err_cdf_ln).^2); %Sum squared error lognormal
            % Get the beta prime approximation CDF
            [k1,k2,theta1,theta2] = beta_prime_approx_v2(hL, hi, dL, di, K, omega, sigma_N, noise_dB, n);
            CDF_sinr_beta = beta_prime_cdf(sinr_th, k1, k2, 1, theta1/theta2);
            CDF_sinr_beta(isnan(CDF_sinr_beta))=0;
            err_cdf_beta = (CDF_sinr_beta-CDF_sinr_sim);
            err_cdf_beta(isnan(err_cdf_beta))=CDF_sinr_beta(isnan(err_cdf_beta));
            err_cdf_beta(isinf(err_cdf_beta))=CDF_sinr_beta(isinf(err_cdf_beta));
            CVM_beta(k,j,i) = sum((err_cdf_beta).^2); %Sum squared error beta prime
            % Get the Gauss Hermite approximation CDF
            Np_load = load("Np_CDF_suburban_1.mat");
            Np1 = round(interp1(Np_load.K_L_dB,Np_load.Np_converge,10*log10(K(1)))); % Num CDF terms
            CDF_sinr_gh = gauss_hermite_rician_cdf4_v2(sinr_th, mu_z_0, sigma_z_0, K(1), omega, Np1);
            CDF_sinr_gh(isnan(CDF_sinr_gh))=0;
            CDF_sinr_gh(isinf(CDF_sinr_gh))=0;
            err_cdf_gh = (CDF_sinr_gh-CDF_sinr_sim);
            err_cdf_gh(isnan(err_cdf_gh))=CDF_sinr_gh(isnan(err_cdf_gh));
            err_cdf_gh(isinf(err_cdf_gh))=CDF_sinr_gh(isinf(err_cdf_gh));
            CVM_gh(k,j,i) = sum((err_cdf_gh).^2); %Sum absolute error Gauss Hermite
        end
    end
end

%%
%Save results to file
save("CVM_Convergence.mat",'CVM_ln','CVM_beta','CVM_gh','K_L_dB','sigma_N_L_dB',...
    'dL_H','di','Pt','Pi','height','freq','noise','theta_GCS','n','num_points');

%%
% Plotting the convergence of CvM
CVM = load("CVM_Convergence.mat");

for i = 1:length(K_L_dB)
    for j = 1:length(sigma_N_L_dB)
        figure();
        plot(num_points, CVM.CVM_ln(i,j,:),'ro','LineWidth',1.2,"MarkerSize",9);
        hold on
        plot(num_points, CVM.CVM_beta(i,j,:),'sb','LineWidth',1.2,"MarkerSize",8);
        hold on
        plot(num_points, CVM.CVM_gh(i,j,:),'gd','LineWidth',1,"MarkerSize",8);
        xlabel("No. Points");
        ylabel("CVM");
        title(sprintf("Convergence of CVM K = %d dB, \sigma = %d dB",K_L_dB(i),sigma_N_L_dB(j)));
    end
end

%%
% Plotting the convergence of delta CvM
CVM = load("CVM_Convergence.mat");

for i = 1:length(K_L_dB)
    for j = 1:length(sigma_N_L_dB)
        delta_CVM_ln = (CVM.CVM_ln(i,j,2:end) - CVM.CVM_ln(i,j,1:end-1)) ./ CVM.CVM_ln(i,j,1:end-1) .*100;
        delta_CVM_beta = (CVM.CVM_beta(i,j,2:end) - CVM.CVM_beta(i,j,1:end-1)) ./ CVM.CVM_beta(i,j,1:end-1) .*100;
        delta_CVM_gh = (CVM.CVM_gh(i,j,2:end) - CVM.CVM_gh(i,j,1:end-1)) ./ CVM.CVM_gh(i,j,1:end-1) .*100;
        figure();
        plot(num_points(2:end), delta_CVM_ln,'ro','LineWidth',1.2,"MarkerSize",9);
        hold on
        plot(num_points(2:end), delta_CVM_beta,'sb','LineWidth',1.2,"MarkerSize",8);
        hold on
        plot(num_points(2:end), delta_CVM_gh,'gd','LineWidth',1,"MarkerSize",8);
        xlabel("No. Points");
        ylabel("\delta CVM");
        title(sprintf("Convergence of \delta CVM K = %d dB, \sigma = %d dB",K_L_dB(i),sigma_N_L_dB(j)));
    end
end

function y = log_normal_cdf(x,mu,sigma)
    y = 0.5 + 0.5.*erf((log(x)-mu)./(sqrt(2)*sigma));
end

function y = beta_prime_cdf(x,a,b,p,q)
    %The CDF of the generalised beta prime of the second kind
    %is the regularized incomplete beta function
    %NOTE: MatLab's betainc is the regularized incomplete beta function
    y = betainc(x.^p./(x.^p+q^p),a,b);
end

function h_pl = friis_cnst(Pt,G_Tx,G_Rx,freq)
    lambda = 3e8/freq;
    h_pl = Pt*(G_Tx*G_Rx*lambda^2)/(16*pi^2);
end