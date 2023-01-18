%Date Modified: 3/12/2021
% U2G Sim with 10^9 samples

%Clear and load the Python func
clear classes;
pyfunc = py.importlib.import_module('gh_rician_cdf_F1');
py.importlib.reload(pyfunc);

% Specify the input parameters ===========================================
Pt_dBm = 40; Pt = 10^(Pt_dBm/10)/1000; %40dBm = 10W Tx power
Pi_dBm = 10; Pi = 10^(Pi_dBm/10)/1000; %23dBm = 0.1995W Tx power
freq = 0.968e9; %Channel frequency
height = 20; %Height btw UAVs and GCS (m)
dL_H = 300; %Horizontal distance between leader UAV and GCS
dL = sqrt(dL_H^2 + height^2);
Num_member_UAV = 3; %Number of member UAVs 
radius = 5; % Distances between each follower UAV and the leader UAV
x_i = zeros(1,Num_member_UAV); % X-coord of each member UAV, with GCS as center
y_i = zeros(1,Num_member_UAV); % Y-coord of each member UAV, with GCS as center
PLoS_i = zeros(1,Num_member_UAV);
theta_i_d = zeros(1,Num_member_UAV);
for i = 1:Num_member_UAV
    x_i(i) = dL_H + radius * cos((i-1)*2*pi/Num_member_UAV);
    y_i(i) = radius * sin((i-1)*2*pi/Num_member_UAV);
    PLoS_i(i) = PLoS_v3(sqrt(x_i(i)^2+y_i(i)^2), height, 0, 0.1 ,7.5e-4, 8);
    theta_i_d(i) = atand(height/sqrt(x_i(i)^2+y_i(i)^2));
end
di = sqrt(x_i.^2 + y_i.^2 + height^2);
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
%Calculate channel parameters affected by positions
PLoS_L = PLoS_v3(dL_H,height,0,0.1,7.5e-4,8); % LoS Prob. btw leader & GCS
theta_L_d = atand(height/dL_H); %Elevation angle of leader from GCS in degrees
% G_GCS_L = antenna_gain(theta_L_d,theta_GCS,Go); % Antenna gain from GCS to leader
% G_L_GCS = antenna_gain(theta_L_d,0,Go); % Antenna gain from leader to GCS
% G_i = antenna_gain(0,0,Go); % Antenna gain between leader and member i
G_GCS_L = 1; G_L_GCS = 1; G_i = 1;
h_GCS = friis_cnst(Pt,G_GCS_L,G_L_GCS,freq); %Channel gain U2G
h_UAV = friis_cnst(Pi,G_i,G_i,freq); %Channel gain U2U
n_L = (n_min-n_max)*PLoS_L + n_max; %Path loss exponent btw leader and GCS
n_i = (n_min-n_max).*PLoS_i + n_max; %Path loss exponent btw member UAV and GCS
sigma_N_L_dB = alpha*exp(-beta*theta_L_d); sigma_N_L = 10^(sigma_N_L_dB/10);%Logarithmic Std dev of shadowing btw leader & GCS
sigma_N_i_dB = alpha.*exp(-beta.*theta_i_d); sigma_N_i = 10.^(sigma_N_i_dB./10);%Logarithmic Std dev of shadowing btw leader & member i
K_L_dB = K_dB_min*exp(log(K_dB_max/K_dB_min)*PLoS_L^2); K_L = 10^(K_L_dB/10); %Rician K factor btw leader and GCS
K_i_dB = K_dB_min.*exp(log(K_dB_max/K_dB_min).*PLoS_i.^2); K_i = 10.^(K_i_dB./10); %Rician K factor btw member and GCS
% K_L_dB = 8; K_L = 10^(K_L_dB/10);
% sigma_N_L_dB = 8; sigma_N_L = 10^(sigma_N_L_dB/10);
%Simulating SINR using Monte Carlo =======================================
N = 10^8; %Number of Monte Carlo data points (divided by 10 of the actual number)
%K = 10^(K_dB/10); v = sqrt(omega*K/(1+K)); s = sqrt(omega/(2+2*K));
v_L = sqrt(omega*K_L/(1+K_L)); s_L = sqrt(omega/(2+2*K_L)); %Rician parameters for leader UAV
v_i = sqrt(omega.*K_i./(1+K_i)); s_i = sqrt(omega./(2+2.*K_i)); %Rician parameters for member UAV
Pwr_sim = h_UAV*dL^(-n_L).*ricePwrrnd(v_L*ones(10, N),s_L).*10.^(normrnd(0,sigma_N_L,[10,N])./10);
I_sim = ones(10, N, Num_member_UAV);
for i = 1:Num_member_UAV
    % Simulating interference power from other member UAVs to GCS
    I_sim(:,:,i) = h_UAV*(di(i)).^(-n_i(i)).*ricePwrrnd(v_i(i)*ones(10, N),s_i(i)).*10.^(normrnd(0,sigma_N_i(i),[10,N])./10);
end
if ~isempty(di)
    I_sim = sum(I_sim,3);
end
SINR_sim = Pwr_sim./(I_sim+noise); %SINR
%Get mean, std dev and histogram of SINR
eta = log(10)/10;
% mean_sinr = mean(SINR_sim); %Mean of SINR
% sigma_sinr = std(SINR_sim); %Std dev of SINR
% pdf_range = mean_sinr + 5*sigma_sinr;
% num_bin = 300; % The number of bins to have between 0 to pdf_range
% bin_width = pdf_range/num_bin;
% h = histogram(SINR_sim,'BinWidth',bin_width,'Normalization','pdf'); hold on;
% % h = histogram(SINR_sim,150,'BinLimits',[max(0,mean_sinr-5*sigma_sinr),mean_sinr+5*sigma_sinr],'Normalization','pdf'); hold on;
% % h = histogram(SINR_sim,150,'Normalization','pdf'); hold on;
% pdf_sim = h.Values(1:num_bin); %PDF values
% c = conv(h.BinEdges, [0.5 0.5], 'valid'); %PDF points
% c = c(1:num_bin);

%Obtaining analytical model using approximations =========================
%Lognormal Approximation
K = [K_L,K_i];
sigma_N = [sigma_N_L,sigma_N_i];
n = [n_L, n_i];
hI = h_UAV * ones(1,Num_member_UAV);
[mu_z_0,sigma_z_0] = log_normal_large_scale_approx_general(h_UAV, hI, dL, di, K(2:end), omega, sigma_N, noise_dB, n);
[mu_z_ln,sigma_z_ln] = lognormal_rice_moment_match(eta*mu_z_0,eta*sigma_z_0,K(1),omega);
% approx_sinr_ln = log_normal_pdf(c, eta*mu_z_ln, eta*sigma_z_ln);
% plot(c, approx_sinr_ln, 'r');

%Beta Prime Approximation
[k1,k2,theta1,theta2] = beta_prime_approx_v2_general(h_UAV, hI, dL, di, K, omega, sigma_N, noise_dB, n);
% approx_sinr_beta = beta_prime_pdf(c, k1, k2, 1, theta1/theta2);
% plot(c, approx_sinr_beta, 'g');

%Gauss-Hermite hybrid method
%Version two uses the outputs of log_normal_large_scale_approx.m
% Np = 32; % Gauss Hermite points for PDF
%The var mu_z_0 is multiplied by -1 below since gauss_hermite_rician_pdf_v2
%negates it again.
% approx_sinr_gh = gauss_hermite_rician_pdf_v2(c, mu_z_0, sigma_z_0, K(1), omega, Np);
% plot(c, approx_sinr_gh, 'k');

% xlabel("S")
% ylabel("PDF")
% title("PDF of SINR")
% %title(sprintf("K = %d dB",K_L_dB));
% % legend("Simulated histogram", "Lognormal","Beta Prime","Gauss-Hermite")
% legend("Simulated", "Lognormal","\beta'")

% Plot CDF of SINR ======================================================
% Load the Np vs K data to determine Np
% Np_load = load("Np_CDF_suburban_1.mat");
% Np = round(interp1(Np_load.K_L_dB,Np_load.Np_converge,10*log10(K(1)))); % Num CDF terms
% h = histogram(SINR_sim,150,'BinLimits',[max(0,mean_sinr-5*sigma_sinr),mean_sinr+5*sigma_sinr],'Normalization','cdf');
% CDF_sinr_sim = h.Values;
% c = conv(h.BinEdges, [0.5 0.5], 'valid'); %CDF points
    
%Get c_index from when CDF_sinr_sim is non-zero
% for i = 1:length(CDF_sinr_sim)
%         if CDF_sinr_sim(i) > 0
%             c_index = i;
%             break;
%         end
% end
% c_points = c(c_index:end); 
% c_points = c(end);

mean_sinr = mean(SINR_sim,'all'); %Mean of SINR
sigma_sinr = std(SINR_sim,0,'all'); %Std dev of SINR
cdf_range = mean_sinr + 5*sigma_sinr;
num_points = 150;
c = linspace(min(SINR_sim,[],'all'),cdf_range,num_points+1);
CDF_sinr_sim = zeros(1,num_points);
for i = 1:num_points
    CDF_sinr_sim(i) = sum(SINR_sim<=c(i+1),'all')/(10*N);
end
c_points = c(2:end);

%Gauss-Hermite Approx of CDF
%CDF_sinr_gh = gauss_hermite_rician_cdf4(c_points, Pt, Pi, dL, di, K, K, sigma_N, noise_dB, n, 32);
%CDF_sinr_gh = gauss_hermite_rician_cdf4_py(c_points, Pt, Pi, dL, di, K, K, sigma_N, noise_dB, n, 32);
% tic
% CDF_sinr_gh = gauss_hermite_rician_cdf4_v2(c_points, mu_z_0, sigma_z_0, K(1), omega, Np);
% avg_time_gh = toc/length(c_points)
tic
CDF_sinr_ln = log_normal_cdf(c_points,eta*mu_z_ln, eta*sigma_z_ln);
avg_time_ln = toc/length(c_points)
tic
CDF_sinr_beta = beta_prime_cdf(c_points,k1, k2, 1, theta1/theta2);
avg_time_beta = toc/length(c_points)

% Plot normal CDF
figure();
plot(c_points,CDF_sinr_sim,'k','LineWidth',2); hold on;
plot(c_points(1:10:end),CDF_sinr_ln(1:10:end),'ro','LineWidth',1.2,"MarkerSize",9); hold on;
plot(c_points(1:10:end),CDF_sinr_beta(1:10:end),'sb','LineWidth',1.2,"MarkerSize",8); hold on;
% plot(c_points,CDF_sinr_gh,'gd','LineWidth',1,"MarkerSize",8);
xlabel("s",'FontSize',36);
ylabel("CDF",'FontSize',36);
title("CDF of SINR",'FontSize',36);
%title(sprintf("K = %d dB",K_L_dB));
legend("Simulated","Lognormal","\beta'",'FontSize',36);
set(gca,"FontSize",28) % To change axis fontsize

max_diff_beta = max(abs(CDF_sinr_beta-CDF_sinr_sim))
max_diff_ln = max(abs(CDF_sinr_ln-CDF_sinr_sim))

%%
mean_beta_prime(k1, k2, 1, theta1/theta2)

function y = log_normal_pdf(x,mu,sigma)
    y = exp(-(log(x)-mu).^2./(2*sigma^2))./(sqrt(2*pi)*sigma.*x);
end

function y = gamma_pdf(x,k,theta)
    y = (x.^(k-1)./(theta^k.*gamma(k))).*exp(-x./theta);
end

function y = beta_prime_pdf(x,a,b,p,q)
    num = p.*(x./q).^(a*p-1).*(1+(x./q).^p).^(-a-b); %Numerator of pdf
    den = q*beta(a,b);
    y = num./den;
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

function m = mean_beta_prime(a,b,p,q)
    m = (q*gamma(b-1/p)*gamma(a+1/p)) / (gamma(a)*gamma(b));
end