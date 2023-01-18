% Date: 21/12/2022
% To test central limit theorem on sum of interferers

Pi_dBm = 10; Pi = 10^(Pi_dBm/10)/1000; %23dBm = 0.1995W Tx power
freq = 0.968e9;
hi = friis_cnst(Pi,1,1,freq); %Channel gain U2U
height = 20; %Height btw UAVs and GCS (m)
dL_H = 200; %Horizontal distance between leader UAV and GCS
Num_member_UAV = 499; %Number of member UAVs 
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
omega = 1;   %Rician Scale factor
n_max = 2.75; n_min = 2; %Range of path loss exponent
K_dB_max = 17.5; K_dB_min = 7.8; %Range of Rician K (dB) in suburban Cleveland
alpha = 11.25; beta = 0.06; %Env parameters for logarithm std dev of shadowing 
n_i = (n_min-n_max).*PLoS_i + n_max; %Path loss exponent btw member UAV and GCS
sigma_N_i_dB = alpha.*exp(-beta.*theta_i_d); sigma_N_i = 10.^(sigma_N_i_dB./10);
K_i_dB = K_dB_min.*exp(log(K_dB_max/K_dB_min).*PLoS_i.^2); K_i = 10.^(K_i_dB./10);
v_i = sqrt(omega.*K_i./(1+K_i)); s_i = sqrt(omega./(2+2.*K_i)); % Rician parameters for member UAV
N = 10^7;
% Generate equavalent Generalised-K RV with same mean and std dev as
% interferers
% gen_K = zeros(Num_member_UAV,N);
L_NCS = zeros(Num_member_UAV,N);
for i = 1:Num_member_UAV
    NCS = ricePwrrnd(v_i(i)*ones(1, N),s_i(i));
    LN = 10.^(normrnd(0,sigma_N_i(i),1,N)./10);
%     mean_NCS = mean(NCS);
%     var_NCS = std(NCS)^2;
%     mean_LN = mean(LN);
%     var_LN = std(LN)^2;
%     k_NCS = mean_NCS^2/var_NCS; % Gamma shape param for NCS RV
%     theta_NCS = var_NCS/mean_NCS; % Gamma scale param for NCS RV
%     k_LN = mean_LN^2/var_LN; % Gamma shape param for NCS RV
%     theta_LN = var_LN/mean_LN; % Gamma scale param for NCS RV
%     gen_K(i,:) = hi*(di(i)).^(-n_i(i)).*gamrnd(k_NCS,theta_NCS,1,N).*gamrnd(k_LN,theta_LN,1,N);
    L_NCS(i,:) = hi*(di(i)).^(-n_i(i)).*NCS.*LN;
end
if Num_member_UAV > 1
%     gen_K = sum(gen_K);
    L_NCS = sum(L_NCS);
end
% Get gamma approx params
m = mean(L_NCS); std_dev = std(L_NCS); var = std_dev^2;
k = m^2/var;
theta = var/m;

% Get the lognormal approx params
% sigma_ln = sqrt(log(1+var/m^2));
% mu_ln = log(m)-sigma_ln^2/2;

figure();
pdf_range = m + 5*std_dev; num_bin = 300; bin_width = pdf_range/num_bin;
% h = histogram(L_NCS,300,'Normalization','pdf'); hold on;
h = histogram(L_NCS,'BinWidth',bin_width,'Normalization','pdf'); hold on;
% h2 = histogram(L_NCS,'BinWidth',bin_width,'Normalization','pdf'); hold on;
NCS_pdf = h.Values; %PDF values
c = conv(h.BinEdges, [0.5 0.5], 'valid'); %PDF points
nomr_pdf = normpdf(c,m,std_dev);
plot(c,nomr_pdf,'r','LineWidth',1)
gamma_PDF = gamma_pdf(c,k,theta,'vpa_mode');
plot(c,gamma_PDF,'g','LineWidth',1.5)
% lognormal_PDF = log_normal_pdf(c,mu_ln,sigma_ln);
% plot(c,lognormal_PDF,'m','LineWidth',1)
xlabel("X"); ylabel("f(X)");
% title(sprintf("K = %d dB",K_dB));
% legend("Simulated","Normal approx");
legend("Simulated","Normal approx","Gamma approx");

figure();
num_points = 150;
c = linspace(min(L_NCS),pdf_range,num_points+1);
NCS_cdf = zeros(1,num_points);
for i = 1:num_points
    NCS_cdf(i) = sum(L_NCS<=c(i+1))/N;
end
c_points = c(2:end);
gamma_CDF = gamma_cdf(c_points,k,theta);
norm_CDF = normcdf(c_points, m, std_dev);
% lognormal_CDF = log_normal_cdf(c_points,mu_ln,sigma_ln);
plot(c_points,NCS_cdf,'LineWidth',1.5); hold on;
plot(c_points(1:2:end),norm_CDF(1:2:end),'or'); hold on;
plot(c_points(1:2:end),gamma_CDF(1:2:end),'sg'); hold on;
% plot(c_points(1:2:end),lognormal_CDF(1:2:end),'om',"MarkerSize",5);
xlabel("x"); ylabel("F(x)");
% legend("Simulated","Normal approx");
legend("Simulated","Normal approx","Gamma approx");

% std_dev_ln = std(log(L_NCS));
% std_dev_ln_dB = 10.*log10(std_dev_ln)
% max_abs_diff = max(abs(gamma_CDF - NCS_cdf))
%%
function h_pl = friis_cnst(Pt,G_Tx,G_Rx,freq)
    lambda = 3e8/freq;
    h_pl = Pt*(G_Tx*G_Rx*lambda^2)/(16*pi^2);
end

function y = gamma_pdf(x,k,theta,vpa_mode)
    if exist('vpa_mode','var')
        x = vpa(x);
        k = vpa(k);
        theta = vpa(theta);
    end
    y = (x.^(k-1)./(theta^k.*gamma(k))).*exp(-x./theta);
end

function y = gamma_cdf(x,k,theta)
    y = gammainc(x./theta,k,'lower');
end