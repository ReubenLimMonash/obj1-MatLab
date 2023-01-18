% Date: 6/1/2023
% Store the CDF of SINR vs height distance using 10e7 sample
% CDF generated using CDF_points from 10e8 samples result
% For G2U communication

%Clear and load the Python func
clear classes;
pyfunc = py.importlib.import_module('gh_rician_cdf_F1');
py.importlib.reload(pyfunc);

% Specify the input parameters ===========================================
Pt_dBm = 40; Pt = 10^(Pt_dBm/10)/1000; %40dBm = 10W Tx power
Pi_dBm = 10; Pi = 10^(Pi_dBm/10)/1000; %23dBm = 0.1995W Tx power
freq = 0.968e9; %Channel frequency
height = 20:20:300; %Height btw UAVs and  GCS (m)
dL_H = 20:20:300; % Horizontal distance between swarm and GCS
Num_member_UAV = 7; %Number of member UAVs
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
N = 10^7; %Number of Monte Carlo data points
Np_load = load("Np_CDF_suburban_1.mat"); %For number of terms in GH method
% num_bin = 150; % The number of bins to have between 0 to cdf_range (if using histogram)
num_points = 300; % The number of points to evaluate CDF at (without using histogram)

% To store the CDF at each combination of height and h dist
CDF_sim = zeros(length(height),length(dL_H),num_points);
CDF_gh = zeros(length(height),length(dL_H),num_points);
CDF_ln = zeros(length(height),length(dL_H),num_points);
CDF_beta = zeros(length(height),length(dL_H),num_points);
CDF_points = zeros(length(height),length(dL_H),num_points);

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
        %G_GCS_L = antenna_gain(theta_L_d,theta_GCS,Go); % Antenna gain from GCS to leader
        %G_L_GCS = antenna_gain(theta_L_d,0,Go); % Antenna gain from leader to GCS
        G_GCS_L = 1; G_L_GCS = 1; G_i = 1;
%         G_i = antenna_gain(0,0,Go); % Antenna gain between leader and member i
        hL = friis_cnst(Pt,G_GCS_L,G_L_GCS,freq); %Channel gain U2G
        hi = friis_cnst(Pi,G_i,G_i,freq); %Channel gain U2U
        n_L = (n_min-n_max)*PLoS + n_max; %Path loss exponent btw leader and GCS
        n_i = (n_min-n_max)*1 + n_max; %Path loss exponent btw leader and member i (PLoS = 1)
        sigma_N_L_dB = alpha*exp(-beta*theta_L_d); sigma_N_L = 10^(sigma_N_L_dB/10);%Logarithmic Std dev of shadowing btw leader & GCS
        sigma_N_i = 10^(sigma_N_i_dB/10);%Logarithmic Std dev of shadowing btw leader & member i
        K_L_dB = K_dB_min*exp(log(K_dB_max/K_dB_min)*PLoS^2); K_L = 10^(K_L_dB/10); %Rician K factor btw leader and GCS
        K_i_dB = K_dB_max; K_i = 10^(K_i_dB/10); %Rician K factor btw leader and member i
        sigma_N_L_dB_store(k,j) = sigma_N_L_dB;
        K_L_dB_store(k,j) = K_L_dB;
        % Monte Carlo sim ====================================================
        v_L = sqrt(omega*K_L/(1+K_L)); s_L = sqrt(omega/(2+2*K_L)); %Rician parameters for leader UAV
        v_i = sqrt(omega*K_i/(1+K_i)); s_i = sqrt(omega/(2+2*K_i)); %Rician parameters for member UAV
        Pwr_sim = hL*dL^(-n_L).*ricePwrrnd(v_L*ones(1, N),s_L).*10.^(normrnd(0,sigma_N_L,1,N)./10);
        I_sim = ones(length(di), N);
        for i = 1:length(di)
            I_sim(i,:) = (di(i)).^(-n_i).*ricePwrrnd(v_i*ones(1, N),s_i).*10.^(normrnd(0,sigma_N_i,1,N)./10);
        end
        if length(di) > 1
            I_sim = hi.*sum(I_sim);
        else
            I_sim = hi.*I_sim;
        end
        SINR_sim = Pwr_sim./(I_sim+noise); %SINR
        %Get simulated SINR CDF
%         mean_sinr = mean(SINR_sim); %Mean of SINR
%         sigma_sinr = std(SINR_sim); %Std dev of SINR
%         cdf_range = mean_sinr + 5*sigma_sinr; 
%         c = linspace(min(SINR_sim),cdf_range,num_points+1);
        c = load("CDF_vs_Height_Dist_10e8/CDF_vs_Height_Dist_7UAV.mat","CDF_points").CDF_points(k,j,:); % Make sure to edit the file name for each no. of UAVs
        CDF_sinr_sim = zeros(1,num_points);
        for i = 1:num_points
%             CDF_sinr_sim(i) = sum(SINR_sim<=c(i+1))/N;
            CDF_sinr_sim(i) = sum(SINR_sim<=c(i))/N;
        end
        sinr_th = c;
        CDF_sim(k,j,:) = CDF_sinr_sim;
        CDF_points(k,j,:) = sinr_th;
        %======================================================================
%         K = [K_L, K_i*ones(1,length(di))]; %Corresponding to [GCS, Member 1, Member 2]
%         sigma_N = [sigma_N_L, sigma_N_i*ones(1,length(di))]; %Corresponding to [GCS, Member 1, Member 2]
%         n = [n_L, n_i*ones(1,length(di))]; %Corresponding to [GCS, Member 1, Member 2]
%         % Get the lognormal approximation CDF
%         [mu_z_0,sigma_z_0] = log_normal_large_scale_approx(hL, hi, dL, di, K(2:end), omega, sigma_N, noise_dB, n);
%         [mu_z_ln,sigma_z_ln] = lognormal_rice_moment_match(eta*mu_z_0,eta*sigma_z_0,K(1),omega);
%         CDF_sinr_ln = log_normal_cdf(sinr_th, eta*mu_z_ln, eta*sigma_z_ln);
%         CDF_sinr_ln(isnan(CDF_sinr_ln))=0;
%         CDF_ln(k,j,:) = CDF_sinr_ln;
%         % Get the beta prime approximation CDF
%         [k1,k2,theta1,theta2] = beta_prime_approx_v2(hL, hi, dL, di, K, omega, sigma_N, noise_dB, n);
%         CDF_sinr_beta = beta_prime_cdf(sinr_th, k1, k2, 1, theta1/theta2);
%         CDF_sinr_beta(isnan(CDF_sinr_beta))=0;
%         CDF_beta(k,j,:) = CDF_sinr_beta;
        
    end
end
sim_state = "Done"
%%
% Save the variable
save("CDF_vs_Height_Dist_7UAV_10e7.mat",'dL_H','di','Pt','Pi','height','Num_member_UAV',...
    'CDF_sim','CDF_points','N');

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