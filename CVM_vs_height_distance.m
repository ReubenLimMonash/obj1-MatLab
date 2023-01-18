% Date: 21/11/2021
% Reliability surface plot vs height and distance

%Clear and load the Python func
clear classes;
pyfunc = py.importlib.import_module('gh_rician_cdf_F1');
py.importlib.reload(pyfunc);

% Specify the input parameters ===========================================
Pt_dBm = 40; Pt = 10^(Pt_dBm/10)/1000; %40dBm = 10W Tx power
Pi_dBm = 10; Pi = 10^(Pi_dBm/10)/1000; %23dBm = 0.1995W Tx power
freq = 0.968e9; %Channel frequency
height = 50:25:300; %Height btw UAVs and  GCS (m)
dL_H = 50:25:300; % Horizontal distance between swarm and GCS
% dL = sqrt(dL_H.^2 + height^2);
Num_member_UAV = 3; %Number of member UAVs
di = 5; %Distances between each follower UAV and the leader UAV
di = ones(1,Num_member_UAV)*di; % [5, 5, 5]
noise_dB = -107;
noise = 10^(noise_dB/10)/1000;
omega = 1;   %Rician Scale factor
theta_GCS = 30; % Antenna tilt angle of GCS (in degrees)
Go = 1; % Antenna gain at horizontal
%Environment parameters-----------------------------
%ENVIRONMENT = SUBURBAN
n_max = 2.75; n_min = 2; %Range of path loss exponent
K_dB_max = 17.5; K_dB_min = 7.8; %Range of Rician K (dB) in suburban Cleveland
alpha = 11.25; beta = 0.06; %Env parameters for logarithm std dev of shadowing
%sigma_N_i_dB = 1.9144; %Logarithmic Std dev of shadowing btw leader & member i (in dB)
sigma_N_i_dB = -inf;
eta = log(10)/10;
N = 10^8; %Number of Monte Carlo data points
Np_load = load("Np_CDF_suburban_1.mat"); %For number of terms in GH method

CVM_ln = zeros(length(height),length(dL_H)); % Cramer-von Mises criterion for LN
CVM_beta = zeros(length(height),length(dL_H)); % Cramer-von Mises criterion for BP
CVM_gh = zeros(length(height),length(dL_H)); % Cramer-von Mises criterion for GH
reliability_sim = zeros(length(height),length(dL_H));
reliability_gh = zeros(length(height),length(dL_H));
reliability_ln = zeros(length(height),length(dL_H));
reliability_beta = zeros(length(height),length(dL_H));

%Let's also store the sigma_N_L and K_L_dB values
sigma_N_L_dB_store = zeros(length(height),length(dL_H));
K_L_dB_store = zeros(length(height),length(dL_H));

for j = 1:length(height)
    height(j)
    for k = 1:length(dL_H)
        dL_H(k)
        dL = sqrt(dL_H(k).^2 + height(j)^2);
        %Calculate channel parameters affected by positions
        PLoS = PLoS_v3(dL_H(k),0,height(j),0.1,7.5e-4,8); % LoS Prob. btw leader & GCS
        theta_L_d = atand(height(j)/dL_H(k)); %Elevation angle of leader from GCS in degrees
        % G_GCS_L = antenna_gain(theta_L_d,theta_GCS,Go); % Antenna gain from GCS to leader
        % G_L_GCS = antenna_gain(theta_L_d,0,Go); % Antenna gain from leader to GCS
        % G_i = antenna_gain(0,0,Go); % Antenna gain between leader and member i
        G_GCS_L = 1; G_L_GCS = 1; G_i = 1; 
        hL = friis_cnst(Pt,G_GCS_L,G_L_GCS,freq); %Channel gain U2G
        hi = friis_cnst(Pi,G_i,G_i,freq); %Channel gain U2U
        n_L = (n_min-n_max)*PLoS + n_max; %Path loss exponent btw leader and GCS
        n_i = (n_min-n_max)*1 + n_max; %Path loss exponent btw leader and member i (PLoS = 1)
        sigma_N_L_dB = alpha*exp(-beta*theta_L_d); sigma_N_L = 10^(sigma_N_L_dB/10);%Logarithmic Std dev of shadowing btw leader & GCS
        sigma_N_i = 10^(sigma_N_i_dB/10);%Logarithmic Std dev of shadowing btw leader & member i
        K_L_dB = K_dB_min*exp(log(K_dB_max/K_dB_min)*PLoS^2); K_L = 10^(K_L_dB/10); %Rician K factor btw leader and GCS
        K_i_dB = K_dB_max; K_i = 10^(K_i_dB/10); %Rician K factor btw leader and member i
        sigma_N_L_dB_store(j,k) = sigma_N_L_dB;
        K_L_dB_store(j,k) = K_L_dB;
        % Monte Carlo sim ====================================================
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
        h = histogram(SINR_sim,150,'BinLimits',[max(0,mean_sinr-5*sigma_sinr),mean_sinr+5*sigma_sinr],'Normalization','cdf');
        % The following h = ... uses TOO LITTLE BINS for the whole SINR range
        %h = histogram(SINR_sim,150,'BinLimits',[max(0,mean_sinr-5*sigma_sinr),max(SINR_sim)],'Normalization','cdf');
        %h = histogram(SINR_sim,150,'BinLimits',[0,max(SINR_sim)],'Normalization','cdf');
        CDF_sinr_sim = h.Values;
        c = conv(h.BinEdges, [0.5 0.5], 'valid'); %CDF points
        c_index = CDF_sinr_sim > 0;
        sinr_th = c(c_index);
        CDF_sinr_sim = CDF_sinr_sim(c_index);
        % reliability_sim(j,k) = 1 - interp1(c,CDF_sinr_sim,sinr_th);
        %======================================================================
        % Approximation of CDF of SINR using Gauss Hermite
        K = [K_L, K_i*ones(1,length(di))]; %Corresponding to [GCS, Member 1, Member 2]
        sigma_N = [sigma_N_L, sigma_N_i*ones(1,length(di))]; %Corresponding to [GCS, Member 1, Member 2]
        n = [n_L, n_i*ones(1,length(di))]; %Corresponding to [GCS, Member 1, Member 2]
        % Get the Gauss Hermite approximation CDF
        Np = round(interp1(Np_load.K_L_dB,Np_load.Np_converge,10*log10(K(1)))); % Num CDF terms
        [mu_z_0,sigma_z_0] = log_normal_large_scale_approx(hL, hi, dL, di, K(2:end), omega, sigma_N, noise_dB, n);
        CDF_sinr_gh = gauss_hermite_rician_cdf4_v2(sinr_th, mu_z_0, sigma_z_0, K(1), omega, Np);
        CDF_sinr_gh(isnan(CDF_sinr_gh))=0;
        CDF_sinr_gh(isinf(CDF_sinr_gh))=0;
        err_cdf_gh = (CDF_sinr_gh-CDF_sinr_sim);
        err_cdf_gh(isnan(err_cdf_gh))=CDF_sinr_gh(isnan(err_cdf_gh));
        err_cdf_gh(isinf(err_cdf_gh))=CDF_sinr_gh(isinf(err_cdf_gh));
        CVM_gh(j,k) = sum((err_cdf_gh).^2); %Sum absolute error Gauss Hermite

        % Approximation of CDF of SINR using lognormal
        [mu_z_ln,sigma_z_ln] = lognormal_rice_moment_match(eta*mu_z_0,eta*sigma_z_0,K(1),omega);
        CDF_sinr_ln = log_normal_cdf(sinr_th,eta*mu_z_ln, eta*sigma_z_ln);
        CDF_sinr_ln(isnan(CDF_sinr_ln))=0;
        err_cdf_ln = (CDF_sinr_ln-CDF_sinr_sim);
        err_cdf_ln(isnan(err_cdf_ln))=CDF_sinr_ln(isnan(err_cdf_ln));
        err_cdf_ln(isinf(err_cdf_ln))=CDF_sinr_ln(isinf(err_cdf_ln));
        CVM_ln(k,j) = sum((err_cdf_ln).^2); %Sum squared error lognormal
        % Approximation of CDF of SINR using beta prime
        [k1,k2,theta1,theta2] = beta_prime_approx_v2(hL, hi, dL, di, K, omega, sigma_N, noise_dB, n);
        CDF_sinr_beta = beta_prime_cdf(sinr_th,k1, k2, 1, theta1/theta2);
        CDF_sinr_beta(isnan(CDF_sinr_beta))=0;
        err_cdf_beta = (CDF_sinr_beta-CDF_sinr_sim);
        err_cdf_beta(isnan(err_cdf_beta))=CDF_sinr_beta(isnan(err_cdf_beta));
        err_cdf_beta(isinf(err_cdf_beta))=CDF_sinr_beta(isinf(err_cdf_beta));
        CVM_beta(j,k) = sum((err_cdf_beta).^2); %Sum squared error beta prime
        %======================================================================
    end
end

%%
% Save the variable
Note = "Antenna gains for U2G set to 1. Pt = 40dBm, Pi = 10dBm. 3 member UAVs";
save("CvM_vs_Height_Dist_Pi10dBm_3UAV.mat",'dL_H','di','Pt','Pi','height','theta_GCS',...
    'CVM_ln','CVM_beta','CVM_gh','K_L_dB_store','sigma_N_L_dB_store','Note');

%%
% Load the saved data
CVM = load("CvM_vs_Height_Dist_Pi10dBm_3UAV.mat");

% Plot using mesh 
X = repmat(CVM.dL_H',1,length(CVM.height));
Y = repmat(CVM.height,length(CVM.dL_H),1);
figure();
m1 = mesh(X,Y,CVM.CVM_gh,'EdgeColor',[0,0.5,0],'FaceColor',[0.25,1,0.25],'FaceAlpha',0.2,...
    'LineStyle','--','LineWidth',1);
% set(m1,'Marker','d');
% set(m1,'MarkerSize',3);
hold on;
m2 = mesh(X,Y,CVM.CVM_beta,'EdgeColor',[0,0,1],'FaceColor',[0,0,1],'FaceAlpha',0.1,...
    'LineStyle','-','LineWidth',1);
% set(m2,'Marker','s');
% set(m2,'MarkerSize',3);
hold on;
m3 = mesh(X,Y,CVM.CVM_ln,'EdgeColor',[1,0,0],'FaceColor',[1,0,0],'FaceAlpha',0.3,...
    'LineStyle',':','LineWidth',1);
% set(m3,'Marker','o');
% set(m3,'MarkerSize',3);
xlabel("Horizontal Distance (m)");
ylabel("Height (m)");
zlabel("\delta");
title("Cramer-von Mises Measure of All Approximations");
legend("L-GH","\beta'","Lognormal");

set(gca, "ZScale", "log")

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