%% 
% Version 2
% Date: 7/10/2022
% Reliability surface plot vs height and distance

%Clear and load the Python func
clear classes;
pyfunc = py.importlib.import_module('gh_rician_cdf_F1');
py.importlib.reload(pyfunc);

% Specify the input parameters ===========================================
Pt_dBm = 40; Pt = 10^(Pt_dBm/10)/1000; %40dBm = 10W Tx power
Pi_dBm = 10; Pi = 10^(Pi_dBm/10)/1000; %23dBm = 0.1995W Tx power
freq = 0.968e9; %Channel frequency
height = 10:10:300; %Height btw UAVs and  GCS (m)
dL_H = 10:10:300; % Horizontal distance between swarm and GCS
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
% M = 32; %Packet size in bytes
B = 5e6; %Bandwidth
% tau = 1e-3; %Transmission delay requirement
% sinr_th = 2^(M*8/(B*tau))-1; % Threshold of SINR
R = 12167;
% R = 1e6;
% R = 5e6;
sinr_th = 2^(R/B) - 1; % Threshold of SINR
Np_load = load("Np_CDF_suburban_1.mat"); %For number of terms in GH method

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
        %G_GCS_L = antenna_gain(theta_L_d,theta_GCS,Go); % Antenna gain from GCS to leader
        %G_L_GCS = antenna_gain(theta_L_d,0,Go); % Antenna gain from leader to GCS
        G_GCS_L = 1; G_L_GCS = 1;
        G_i = antenna_gain(0,0,Go); % Antenna gain between leader and member i
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
        if length(di) > 1
            I_sim = hi.*sum(I_sim);
        else
            I_sim = hi.*I_sim;
        end
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
        % Check if the SINR threshold is above max SINR_sim, if it is the
        % threshold is beyond the range of the CDF, and we don't have data
        % to compute reliability at this point. 
        if sinr_th > max(c)
            % Set the reliability for the rest of the horizontal points to
            % NaN. Break out of evaluating further horizontal distance
            reliability_sim(j,k:end) = NaN;
            reliability_gh(j,k:end) = NaN;
            reliability_ln(j,k:end) = NaN;
            reliability_beta(j,k:end) = NaN;
            break
        end
        reliability_sim(j,k) = 1 - interp1(c,CDF_sinr_sim,sinr_th,'spline'); % Use spline method for interpolation. Extrapolation also supported
        %======================================================================
        % Approximation of CDF of SINR using Gauss Hermite
        K = [K_L, K_i*ones(1,length(di))]; %Corresponding to [GCS, Member 1, Member 2]
        sigma_N = [sigma_N_L, sigma_N_i*ones(1,length(di))]; %Corresponding to [GCS, Member 1, Member 2]
        n = [n_L, n_i*ones(1,length(di))]; %Corresponding to [GCS, Member 1, Member 2]
        % Get the Gauss Hermite approximation CDF
        Np = round(interp1(Np_load.K_L_dB,Np_load.Np_converge,10*log10(K(1)))); % Num CDF terms
        [mu_z_0,sigma_z_0] = log_normal_large_scale_approx(hL, hi, dL, di, K(2:end), omega, sigma_N, noise_dB, n);
        CDF_sinr_gh = gauss_hermite_rician_cdf4_v2(sinr_th, mu_z_0, sigma_z_0, K(1), omega, Np);
        reliability_gh(j,k) = 1 - CDF_sinr_gh;

        % Approximation of CDF of SINR using lognormal
        [mu_z_ln,sigma_z_ln] = lognormal_rice_moment_match(eta*mu_z_0,eta*sigma_z_0,K(1),omega);
        CDF_sinr_ln = log_normal_cdf(sinr_th,eta*mu_z_ln, eta*sigma_z_ln);
        reliability_ln(j,k) = 1 - CDF_sinr_ln;

        % Approximation of CDF of SINR using lognormal
        [k1,k2,theta1,theta2] = beta_prime_approx_v2(hL, hi, dL, di, K, omega, sigma_N, noise_dB, n);
        CDF_sinr_beta = beta_prime_cdf(sinr_th,k1, k2, 1, theta1/theta2);
        reliability_beta(j,k) = 1 - CDF_sinr_beta;
        %======================================================================
    end
end
%Post-processing and calculate errors
% reliability_sim(isnan(reliability_sim)) = 0;
% assert(sum(sum(isnan(reliability_sim)))<=0,"NaN exist in reliability_sim");
% reliability_sim(reliability_sim > 1) = 1;
% reliability_sim(reliability_sim < 0) = 0;
% reliability_sim=zeros(11,11);
abs_err_ln = abs(reliability_ln-reliability_sim);
abs_err_beta = abs(reliability_beta-reliability_sim);
abs_err_gh = abs(reliability_gh-reliability_sim);
rel_err_ln = abs(reliability_ln-reliability_sim)./reliability_sim;
rel_err_beta = abs(reliability_beta-reliability_sim)./reliability_sim;
rel_err_gh = abs(reliability_gh-reliability_sim)./reliability_sim;
%%
figure()
plot(dL_H,reliability_sim);
hold on
plot(dL_H,reliability_ln,'ro');
hold on
plot(dL_H,reliability_beta,'sg');
hold on
plot(dL_H,reliability_gh,'dk');
xlabel("Distance from GCS (m)");
ylabel("Reliability");
%title("Reliability of UAV communication")
legend("Simulated","Lognormal","Beta prime","Gauss-Hermite");
%%
% Save the variable
Note = "Antenna gains for U2G set to 1. R = 12.167kbps. 3 member UAV";
save("Reliability_vs_Height_Dist_Pi10dBm_R12167_3UAV.mat",'dL_H','di','Pt','Pi','height','theta_GCS',...
    'reliability_ln','reliability_beta','reliability_gh','reliability_sim',...
    'abs_err_ln','abs_err_beta','abs_err_gh','rel_err_ln','rel_err_beta','rel_err_gh','K_L_dB_store','sigma_N_L_dB_store','Note');

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