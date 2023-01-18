 %Date: 18/10/2021
% To calculate number of terms needed in Gauss-Hermite method for CDF
% for different values of K

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
noise_dBm = -107;
noise = 10^(noise_dBm/10)/1000;
omega = 1;   %Rician Scale factor
theta_GCS = 45; % Antenna tilt angle of GCS (in degrees)
Go = 1; % Antenna gain at horizontal
%Environment parameters-----------------------------
%ENVIRONMENT = SUBURBAN
n_max = 2.75; n_min = 2; %Range of path loss exponent
K_dB_max = 17.5; K_dB_min = 7.8; %Range of Rician K (dB) in suburban Cleveland
alpha = 11.25; beta = 0.06; %Env parameters for logarithm std dev of shadowing 
%sigma_N_i_dB = 1.9144; %Logarithmic Std dev of shadowing btw leader & member i (in dB)
sigma_N_i_dB = -inf;
%Calculate channel parameters affected by positions
PLoS = PLoS_v3(dL_H,height,0,0.1,7.5e-4,8); % LoS Prob. btw leader & GCS
theta_L_d = atand(height/dL_H); %Elevation angle of leader from GCS in degrees
G_GCS_L = antenna_gain(theta_L_d,theta_GCS,Go); % Antenna gain from GCS to leader
G_L_GCS = antenna_gain(theta_L_d,0,Go); % Antenna gain from leader to GCS
G_i = antenna_gain(0,0,Go); % Antenna gain between leader and member i
hL = friis_cnst(Pt,G_GCS_L,G_L_GCS,freq); %Channel gain U2G
hi = friis_cnst(Pi,G_i,G_i,freq); %Channel gain U2U
n_L = (n_min-n_max)*PLoS + n_max; %Path loss exponent btw leader and GCS
n_i = (n_min-n_max)*1 + n_max; %Path loss exponent btw leader and member i (PLoS = 1)
sigma_N_L_dB = alpha*exp(-beta*theta_L_d); sigma_N_L = 10^(sigma_N_L_dB/10);%Logarithmic Std dev of shadowing btw leader & GCS
sigma_N_i = 10^(sigma_N_i_dB/10);%Logarithmic Std dev of shadowing btw leader & member i
%K_L_dB = K_dB_min*exp(log(K_dB_max/K_dB_min)*PLoS^2); K_L = 10^(K_L_dB/10); %Rician K factor btw leader and GCS
K_i_dB = K_dB_max; K_i = 10^(K_i_dB/10); %Rician K factor btw leader and member i
eta = log(10)/10;
N = 10^8; %Number of Monte Carlo data points

K_L_dB = 8:1:17;
K_L_dB = [7.8, K_L_dB, 17.5];
%K_L_dB = 18;
%K_L_dB = [17, 17.5];
Np_converge = ones(1,length(K_L_dB)); % Value of Np where GH converges, for each K
Np_start = 10; %The start value of Np to use (will be updated with prev. result)
for k = 1:length(K_L_dB)
    %sigma_N(k)
    % First, lets generate the simulated values using Monte Carlo sim
    K_L = 10^(K_L_dB(k)/10); K_L_dB(k)
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
    CDF_sinr_sim = h.Values;
    c = conv(h.BinEdges, [0.5 0.5], 'valid'); %CDF points
    %Find the first index where CDF_sinr_sim is > 0
    for i = 1:length(CDF_sinr_sim)
        if CDF_sinr_sim(i) > 0
            c_index = i;
            break;
        end
    end
    %======================================================================
    % Approximation of CDF of SINR using Gauss Hermite
    K = [K_L, K_i, K_i, K_i]; %Corresponding to [GCS, Member 1, Member 2]
    sigma_N = [sigma_N_L, sigma_N_i, sigma_N_i, sigma_N_i]; %Corresponding to [GCS, Member 1, Member 2]
    n = [n_L, n_i, n_i, n_i]; %Corresponding to [GCS, Member 1, Member 2]
    % Get the Gauss Hermite approximation CDF
    PE_cdf_gh = [];
    for Np = Np_start:150
        [mu_z_0,sigma_z_0] = log_normal_large_scale_approx(hL, hi, dL, di, K(2:end), omega, sigma_N, noise_dBm, n);
        CDF_sinr_gh = gauss_hermite_rician_cdf4_v2(c(c_index), mu_z_0, sigma_z_0, K(1), omega, Np);
        if isnan(CDF_sinr_gh) || isinf(CDF_sinr_gh)
            CDF_sinr_gh = 1; %If this cond true, Np is too low, hence default CDF to 1
                             %Error will be high since we are taking the
                             %since we are taking the lowest CDF_sim value
        end
        PE_cdf_gh(end+1) = abs(CDF_sinr_gh-CDF_sinr_sim(c_index))./CDF_sinr_sim(c_index)*100; %Percentage error Gauss Hermite
        % If the slope of the percentage error graph is within 10,
        % take that point as convergence
        if length(PE_cdf_gh) > 1
            delta_PE = PE_cdf_gh(end) - PE_cdf_gh(end-1);
            if (abs(delta_PE) < 5)
                Np_converge(k) = Np; 
                Np_start = Np; %Update Np_start to narrow down search of next Np
                                 %Increasing K needs Np to increase, so
                                 %next iter start from latest Np_converge
                Np_converge(k)   %Display the value for progress
                break;
            end
        end
    end
    %======================================================================
end
%%
figure()
plot(K_L_dB,Np_converge);
xlabel("K (dB)");
ylabel("Number of terms");
title("Num terms for different K")

% Save the variable
% K_dB = K_L_dB;
% state = 'suburban';
% save("Np_CDF_suburban_1.mat",'Np_converge','K_L_dB','dL_H','di','height',...
%     'sigma_N_L_dB','noise_dBm','Pi','Pt','freq','state');

%%
%For Timing Measurements
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
noise_dBm = -107;
noise = 10^(noise_dBm/10)/1000;
omega = 1;   %Rician Scale factor
theta_GCS = 45; % Antenna tilt angle of GCS (in degrees)
Go = 1; % Antenna gain at horizontal
%Environment parameters-----------------------------
%ENVIRONMENT = SUBURBAN
n_max = 2.75; n_min = 2; %Range of path loss exponent
K_dB_max = 17.5; K_dB_min = 7.8; %Range of Rician K (dB) in suburban Cleveland
alpha = 11.25; beta = 0.06; %Env parameters for logarithm std dev of shadowing 
%sigma_N_i_dB = 1.9144; %Logarithmic Std dev of shadowing btw leader & member i (in dB)
sigma_N_i_dB = -inf;
%Calculate channel parameters affected by positions
PLoS = PLoS_v3(dL_H,height,0,0.1,7.5e-4,8); % LoS Prob. btw leader & GCS
theta_L_d = atand(height/dL_H); %Elevation angle of leader from GCS in degrees
G_GCS_L = antenna_gain(theta_L_d,theta_GCS,Go); % Antenna gain from GCS to leader
G_L_GCS = antenna_gain(theta_L_d,0,Go); % Antenna gain from leader to GCS
G_i = antenna_gain(0,0,Go); % Antenna gain between leader and member i
hL = friis_cnst(Pt,G_GCS_L,G_L_GCS,freq); %Channel gain U2G
hi = friis_cnst(Pi,G_i,G_i,freq); %Channel gain U2U
n_L = (n_min-n_max)*PLoS + n_max; %Path loss exponent btw leader and GCS
n_i = (n_min-n_max)*1 + n_max; %Path loss exponent btw leader and member i (PLoS = 1)
sigma_N_L_dB = alpha*exp(-beta*theta_L_d); sigma_N_L = 10^(sigma_N_L_dB/10);%Logarithmic Std dev of shadowing btw leader & GCS
sigma_N_i = 10^(sigma_N_i_dB/10);%Logarithmic Std dev of shadowing btw leader & member i
%K_L_dB = K_dB_min*exp(log(K_dB_max/K_dB_min)*PLoS^2); K_L = 10^(K_L_dB/10); %Rician K factor btw leader and GCS
K_i_dB = K_dB_max; K_i = 10^(K_i_dB/10); %Rician K factor btw leader and member i
eta = log(10)/10;
N = 10^8; %Number of Monte Carlo data points

load('Np_CDF_suburban_1.mat'); %Load the recorded values
time_rec = zeros(1,length(K_L_dB));
for k = 1:length(K_L_dB)
    %sigma_N(k)
    % First, lets generate the simulated values using Monte Carlo sim
    K_L = 10^(K_L_dB(k)/10); K_L_dB(k)
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
    CDF_sinr_sim = h.Values;
    c = conv(h.BinEdges, [0.5 0.5], 'valid'); %CDF points
    %Find the first index where CDF_sinr_sim is > 0
    for i = 1:length(CDF_sinr_sim)
        if CDF_sinr_sim(i) > 0
            c_index = i;
            break;
        end
    end
    %======================================================================
    % Approximation of CDF of SINR using Gauss Hermite
    K = [K_L, K_i, K_i, K_i]; %Corresponding to [GCS, Member 1, Member 2]
    sigma_N = [sigma_N_L, sigma_N_i, sigma_N_i, sigma_N_i]; %Corresponding to [GCS, Member 1, Member 2]
    n = [n_L, n_i, n_i, n_i]; %Corresponding to [GCS, Member 1, Member 2]
    % Get the Gauss Hermite approximation CDF
    % Run once as a warm up
    [mu_z_0,sigma_z_0] = log_normal_large_scale_approx(hL, hi, dL, di, K(2:end), omega, sigma_N, noise_dBm, n);
    CDF_sinr_gh = gauss_hermite_rician_cdf4_v2(c(c_index), mu_z_0, sigma_z_0, K(1), omega, Np_converge(k));
    tic
    [mu_z_0,sigma_z_0] = log_normal_large_scale_approx(hL, hi, dL, di, K(2:end), omega, sigma_N, noise_dBm, n);
    CDF_sinr_gh = gauss_hermite_rician_cdf4_v2(c(c_index), mu_z_0, sigma_z_0, K(1), omega, Np_converge(k));
    time_rec(k) = toc; % Record the time needed for L-GH method with Np_converge(k) num terms
    %======================================================================
end

% save("Np_CDF_suburban_1.mat",'Np_converge','K_L_dB','dL_H','di','height',...
%     'sigma_N_L_dB','noise_dBm','Pi','Pt','freq','state','time_rec');


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