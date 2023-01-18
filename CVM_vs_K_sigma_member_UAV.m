%Date: 17/10/2022
% Calculating MAE and MRE for CDF approximations of SINR
% w.r.t. K and sigma in 2D plot for U2U comm.
% Leader UAV is receiving from the member UAV

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
dL = sqrt(dL_H^2 + height^2); % Distance between leader UAV and GCS
di = [5 5 5]; %Distances between each follower UAV and the leader UAV (for interfering UAVs only). Two values means 3 member UAVs
dm = 5; % Distance between the member UAV of interest and the leader UAV
noise_dB = -107;
noise = 10^(noise_dB/10)/1000;
omega = 1;   %Rician Scale factor
% theta_GCS = 45; % Antenna tilt angle of GCS (in degrees)
Go = 1; % Antenna gain at horizontal
%Environment parameters-----------------------------
%ENVIRONMENT = SUBURBAN
n_max = 2.75; n_min = 2; %Range of path loss exponent
K_dB_max = 17.5; K_dB_min = 7.8; %Range of Rician K (dB) in suburban Cleveland
alpha = 11.25; beta = 0.06; %Env parameters for logarithm std dev of shadowing 
%sigma_N_i_dB = 1.9144; %Logarithmic Std dev of shadowing btw leader & member i (in dB)
sigma_N_i_dB = -inf;
%Calculate channel parameters affected by positions
PLoS = PLoS_v3(dL_H,0,height,0.1,7.5e-4,8); % LoS Prob. btw leader & GCS
theta_L_d = atand(height/dL_H); %Elevation angle of leader from GCS in degrees
% G_GCS_L = antenna_gain(theta_L_d,theta_GCS,Go); % Antenna gain from GCS to leader
% G_L_GCS = antenna_gain(theta_L_d,0,Go); % Antenna gain from leader to GCS
% G_i = antenna_gain(0,0,Go); % Antenna gain between leader and member i
G_GCS_L = 1; G_L_GCS = 1; G_i = 1;
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
% sigma_N_L_dB = 0; sigma_N_L = 10^(sigma_N_L_dB/10);
%K_L_dB = 14; K_L = 10^(K_L_dB/10); 

K_L_dB = 8:1:16;
% K_L_dB = [7.8, K_L_dB];
sigma_N_L_dB = 0:8;
% sigma_N_L_dB = 0:11;
% sigma_N_L_dB = [sigma_N_L_dB, 11.25];
CVM_ln = zeros(length(K_L_dB),length(sigma_N_L_dB)); % Cramer-von Mises criterion for LN
CVM_beta = zeros(length(K_L_dB),length(sigma_N_L_dB)); % Cramer-von Mises criterion for BP
CVM_gh = zeros(length(K_L_dB),length(sigma_N_L_dB)); % Cramer-von Mises criterion for GH
MAE_ln = zeros(length(K_L_dB),length(sigma_N_L_dB)); % Mean abs error for LN
MAE_beta = zeros(length(K_L_dB),length(sigma_N_L_dB)); % Mean abs error for BP
MAE_gh = zeros(length(K_L_dB),length(sigma_N_L_dB)); % Mean abs error for GH
MRE_ln = zeros(length(K_L_dB),length(sigma_N_L_dB)); % Mean rel error for LN
MRE_beta = zeros(length(K_L_dB),length(sigma_N_L_dB)); % Mean rel error for BP
MRE_gh = zeros(length(K_L_dB),length(sigma_N_L_dB)); % Mean rel error for GH

for k = 1:length(K_L_dB)
    for j = 1:length(sigma_N_L_dB)
        K_L = 10^(K_L_dB(k)/10); 
        sigma_N_L = 10^(sigma_N_L_dB(j)/10); 
        sprintf("K = %d dB, \x3C3 = %d dB",K_L_dB(k),sigma_N_L_dB(j));
        v_L = sqrt(omega*K_L/(1+K_L)); s_L = sqrt(omega/(2+2*K_L)); %Rician parameters for leader UAV
        v_i = sqrt(omega*K_i/(1+K_i)); s_i = sqrt(omega/(2+2*K_i)); %Rician parameters for member UAV
        Pwr_sim = hi*dm^(-n_i).*ricePwrrnd(v_i*ones(1, N),s_i).*10.^(normrnd(0,sigma_N_i,1,N)./10);
        I_sim = ones(length(di)+1, N);
        if ~isempty(di)
            for i = 1:length(di)
                % Simulating interference power from other member UAVs to
                % leader UAV
                I_sim(i,:) = hi*(di(i)).^(-n_i).*ricePwrrnd(v_i*ones(1, N),s_i).*10.^(normrnd(0,sigma_N_i,1,N)./10);
            end
        end
        % Simulating interference power from GCS to leader UAV
        I_sim(length(di)+1,:) = hL*dL.^(-n_L).*ricePwrrnd(v_L*ones(1, N),s_L).*10.^(normrnd(0,sigma_N_L,1,N)./10);
        if ~isempty(di)
            I_sim = sum(I_sim);
        end
        SINR_sim = Pwr_sim./(I_sim+noise); %SINR
        %Get simulated SINR CDF
        mean_sinr = mean(SINR_sim); %Mean of SINR
        sigma_sinr = std(SINR_sim); %Std dev of SINR
        h = histogram(SINR_sim,150,'BinLimits',[max(0,mean_sinr-5*sigma_sinr),mean_sinr+5*sigma_sinr],'Normalization','cdf');
%         h = histogram(SINR_sim,300,'BinLimits',[0,max(SINR_sim)],'Normalization','cdf');
        CDF_sinr_sim = h.Values;
        c = conv(h.BinEdges, [0.5 0.5], 'valid'); %CDF points
        c_index = CDF_sinr_sim > 0;
        sinr_th = c(c_index); 
        CDF_sinr_sim = CDF_sinr_sim(c_index);
        %===============================================_member_UAV=======================
        % Approximations of SINR
        K = [K_i, K_i*ones(1,length(di)), K_L]; %Corresponding to [GCS, Member 1, Member 2]
        sigma_N = [sigma_N_i, sigma_N_i*ones(1,length(di)), sigma_N_L]; %Corresponding to [GCS, Member 1, Member 2]
        n = [n_i, n_i*ones(1,length(di)), n_L]; %Corresponding to [GCS, Member 1, Member 2]
        hI = [hi*ones(1,length(di)), hL];
        % Get the lognormal approximation CDF
        [mu_z_0,sigma_z_0] = log_normal_large_scale_approx_general(hi, hI, dm, [di, dL], K(2:end), omega, sigma_N, noise_dB, n);
        [mu_z_ln,sigma_z_ln] = lognormal_rice_moment_match(eta*mu_z_0,eta*sigma_z_0,K(1),omega);
        CDF_sinr_ln = log_normal_cdf(sinr_th, eta*mu_z_ln, eta*sigma_z_ln);
        CDF_sinr_ln(isnan(CDF_sinr_ln))=0;
        err_cdf_ln = (CDF_sinr_ln-CDF_sinr_sim);
        err_cdf_ln(isnan(err_cdf_ln))=CDF_sinr_ln(isnan(err_cdf_ln));
        err_cdf_ln(isinf(err_cdf_ln))=CDF_sinr_ln(isinf(err_cdf_ln));
        CVM_ln(k,j) = sum((err_cdf_ln).^2); %Sum squared error lognormal
        MAE_ln(k,j) = mean(abs(err_cdf_ln)); %Mean abs error lognormal
        MRE_ln(k,j) = mean(abs(err_cdf_ln)./CDF_sinr_sim); %Mean relative error lognormal
        % Get the beta prime approximation CDF
        [k1,k2,theta1,theta2] = beta_prime_approx_v2_general(hi, hI, dm, [di, dL], K, omega, sigma_N, noise_dB, n);
        CDF_sinr_beta = beta_prime_cdf(sinr_th, k1, k2, 1, theta1/theta2);
        CDF_sinr_beta(isnan(CDF_sinr_beta))=0;
        err_cdf_beta = (CDF_sinr_beta-CDF_sinr_sim);
        err_cdf_beta(isnan(err_cdf_beta))=CDF_sinr_beta(isnan(err_cdf_beta));
        err_cdf_beta(isinf(err_cdf_beta))=CDF_sinr_beta(isinf(err_cdf_beta));
        CVM_beta(k,j) = sum((err_cdf_beta).^2); %Sum squared error beta prime
        MAE_beta(k,j) = mean(abs(err_cdf_beta)); %Sum abs error beta prime
        MRE_beta(k,j) = mean(abs(err_cdf_beta)./CDF_sinr_sim); %Mean relative error beta prime
        % Get the Gauss Hermite approximation CDF
%         Np_load = load("Np_CDF_suburban_1.mat");
%         Np1 = round(interp1(Np_load.K_L_dB,Np_load.Np_converge,10*log10(K(1)))); % Num CDF terms
%         CDF_sinr_gh = gauss_hermite_rician_cdf4_v2(sinr_th, mu_z_0, sigma_z_0, K(1), omega, Np1);
%         CDF_sinr_gh(isnan(CDF_sinr_gh))=0;
%         CDF_sinr_gh(isinf(CDF_sinr_gh))=0;
%         err_cdf_gh = (CDF_sinr_gh-CDF_sinr_sim);
%         err_cdf_gh(isnan(err_cdf_gh))=CDF_sinr_gh(isnan(err_cdf_gh));
%         err_cdf_gh(isinf(err_cdf_gh))=CDF_sinr_gh(isinf(err_cdf_gh));
%         CVM_gh(k,j) = sum((err_cdf_gh).^2); %Sum absolute error Gauss Hermite
%         MAE_gh(k,j) = mean(abs(err_cdf_gh)); %Sum abs error Gauss Hermite
%         MRE_gh(k,j) = mean(abs(err_cdf_gh)./CDF_sinr_sim); %Mean relative error Gauss Hermite
    end
end

%%
%Save results to file (for sum of squared error)
save("CVM_vs_K_sigma_1_member_UAV.mat",'CVM_ln','CVM_beta','K_L_dB','sigma_N_L_dB',...
    'dL_H','di','dm','Pt','Pi','height','freq','noise','n');
%%
%Save results to file (for mean of abs error)
save("MAE_vs_K_sigma_1_member_UAV.mat",'MAE_ln','MAE_beta','K_L_dB','sigma_N_L_dB',...
    'dL_H','di','dm','Pt','Pi','height','freq','noise','n');
%%
%Save results to file (for mean of relative error)
save("MRE_vs_K_sigma_1_member_UAV.mat",'MRE_ln','MRE_beta','K_L_dB','sigma_N_L_dB',...
    'dL_H','di','dm','Pt','Pi','height','freq','noise','n');
%%
% Load and plot the results
CVM = load("CVM_vs_K_sigma.mat");
X = repmat(CVM.K_L_dB',1,length(CVM.sigma_N_L_dB));
Y = repmat(CVM.sigma_N_L_dB,length(CVM.K_L_dB),1);
figure();
surf(X,Y,CVM.CVM_gh);
xlabel("K_1_G (dB)");
ylabel("\sigma_1_G (dB)");
zlabel("\delta");
title("Cramer-von Mises Measure for L-GH Approximation");
colorbar('southoutside')

figure();
surf(X,Y,CVM.CVM_beta);
xlabel("K_1_G (dB)");
ylabel("\sigma_1_G (dB)");
zlabel("\delta");
title("Cramer-von Mises Measure for \beta' Approximation");
colorbar('southoutside')

figure();
surf(X,Y,CVM.CVM_ln);
xlabel("K_1_G (dB)");
ylabel("\sigma_1_G (dB)");
zlabel("\delta");
title("Cramer-von Mises Measure for Lognormal Approximation");
colorbar('southoutside')
%%
% Plot using mesh with colour
CVM = load("CVM_vs_K_sigma.mat");
% Plot using mesh 
X = repmat(CVM.K_L_dB',1,length(CVM.sigma_N_L_dB));
Y = repmat(CVM.sigma_N_L_dB,length(CVM.K_L_dB),1);
figure();
m2 = mesh(X,Y,CVM.CVM_beta,'EdgeColor',[0,0,1],'FaceColor',[0,0,1],'FaceAlpha',0.1,...
    'LineStyle','-','LineWidth',1);
% set(m2,'Marker','s');
% set(m2,'MarkerSize',3);
hold on;
m1 = mesh(X,Y,CVM.CVM_gh,'EdgeColor',[0,0.5,0],'FaceColor',[0.25,1,0.25],'FaceAlpha',0.2,...
    'LineStyle','--','LineWidth',1);
% set(m1,'Marker','d');
% set(m1,'MarkerSize',3);
hold on;
m3 = mesh(X,Y,CVM.CVM_ln,'EdgeColor',[1,0,0],'FaceColor',[1,0,0],'FaceAlpha',0.3,...
    'LineStyle',':','LineWidth',1);
% set(m3,'Marker','o');
% set(m3,'MarkerSize',3);
xlabel("K_1_G (dB)");
ylabel("\sigma_1_G (dB)");
zlabel("\delta");
title("Cramer-von Mises Measure of All Approximations");
legend("\beta'","L-GH","Lognormal");



%%
% Plot using mesh with colour
CVM = load("CVM_vs_K_sigma.mat");
X = repmat(CVM.K_L_dB',1,length(CVM.sigma_N_L_dB));
Y = repmat(CVM.sigma_N_L_dB,length(CVM.K_L_dB),1);
figure();
m1 = mesh(X,Y,CVM.CVM_gh,'EdgeColor','g','FaceColor','g','FaceAlpha',0.2,...
    'LineStyle',':','LineWidth',1.3);
set(m1,'Marker','d');
set(m1,'MarkerSize',3);
hold on;
m2 = mesh(X,Y,CVM.CVM_beta,'EdgeColor','b','FaceColor','b','FaceAlpha',0.2,...
    'LineStyle','-','LineWidth',1);
set(m2,'Marker','s');
set(m2,'MarkerSize',3);
hold on;
m3 = mesh(X,Y,CVM.CVM_ln,'EdgeColor','r','FaceColor','r','FaceAlpha',0.2,...
    'LineStyle','--','LineWidth',1.5);
set(m3,'Marker','o');
set(m3,'MarkerSize',3);
hold on;
% plot3([9,9], [1,1], [0,0.4], 'k','LineWidth',1.5);
% hold on;
% plot3([9,9], [7,7], [0,0.4], 'k','LineWidth',1.5);
xlabel("K_1_G (dB)");
ylabel("\sigma_1_G (dB)");
zlabel("\delta");
title("Cramer-von Mises Measure of All Approximations");
legend("L-GH","\beta'","Lognormal");
%%
% Plot using mesh in grayscale
CVM = load("CVM_vs_K_sigma.mat");
X = repmat(CVM.K_L_dB',1,length(CVM.sigma_N_L_dB));
Y = repmat(CVM.sigma_N_L_dB,length(CVM.K_L_dB),1);
figure();
m1 = mesh(X,Y,CVM.CVM_gh,'EdgeColor','k','FaceColor','k','FaceAlpha',0,...
    'LineStyle','-','LineWidth',1);
set(m1,'Marker','d');
set(m1,'MarkerSize',3);
hold on;
m2 = mesh(X,Y,CVM.CVM_beta,'EdgeColor',[0.5,0.5,0.5],'FaceColor','g','FaceAlpha',0,...
    'LineStyle',':','LineWidth',1.3);
set(m2,'Marker','s');
set(m2,'MarkerSize',3);
hold on;
m3 = mesh(X,Y,CVM.CVM_ln,'EdgeColor',[0.25,0.25,0.25],'FaceColor','r','FaceAlpha',0,...
    'LineStyle','--','LineWidth',1.5);
set(m3,'Marker','o');
set(m3,'MarkerSize',3);
hold on;
plot3([9,9], [1,1], [0,0.4], 'k','LineWidth',1.5);
hold on;
plot3([9,9], [7,7], [0,0.4], 'k','LineWidth',1.5);
xlabel("K_1_G (dB)");
ylabel("\sigma_1_G (dB)");
zlabel("\delta");
title("Cramer-von Mises Measure of All Approximations");
legend("L-GH","\beta'","Lognormal");
%%

figure()
plot(sigma_N_L_dB(1:8),smoothdata(CVM_ln(1:8),'sgolay',5),sigma_N_L_dB(1:8),smoothdata(CVM_beta(1:8),'sgolay',5),sigma_N_L_dB(1:8),smoothdata(CVM_gh(1:8),'sgolay',5),'g');
%plot(K_L_dB,MAE_cdf_ln,K_L_dB,MAE_cdf_beta,K_L_dB,MAE_cdf_gh);
hold on
plot(sigma_N_L_dB(1:8),CVM_ln(1:8),'bo','MarkerSize',5);
hold on
plot(sigma_N_L_dB(1:8),CVM_beta(1:8),'ro','MarkerSize',5);
hold on
plot(sigma_N_L_dB(1:8),CVM_gh(1:8),'go','MarkerSize',5);
xlabel("\sigma_1_G (dB)");
ylabel("SAE");
title("Sum of Absolute Error of CDF");
legend("","","","Lognormal","Beta Prime", "Gauss Hermite");

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