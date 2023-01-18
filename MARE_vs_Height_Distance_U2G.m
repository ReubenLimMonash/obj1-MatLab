% Date: 8/11/2022
% Mean Abs Err and Mean Rel Err vs height and distance
% For the communication from the gateway UAV to the GCS

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
% dL = sqrt(dL_H.^2 + height^2);
Num_member_UAV = 31; %Number of member UAVs
radius = 5; % Distances between each follower UAV and the leader UAV
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
% sigma_N_i_dB = -inf;
eta = log(10)/10;
N = 10^8; %Number of Monte Carlo data points
Np_load = load("Np_CDF_suburban_1.mat"); %For number of terms in GH method

% Friis path loss param
%G_GCS_L = antenna_gain(theta_L_d,theta_GCS,Go); % Antenna gain from GCS to leader
%G_L_GCS = antenna_gain(theta_L_d,0,Go); % Antenna gain from leader to GCS
G_GCS_L = 1; G_L_GCS = 1; G_i = 1;
% G_i = antenna_gain(0,0,Go); % Antenna gain between leader and member i
h_GCS = friis_cnst(Pt,G_GCS_L,G_L_GCS,freq); %Channel gain GCS transmitter
h_UAV = friis_cnst(Pi,G_i,G_i,freq); %Channel gain UAV transmitter

MAE_gh = zeros(length(height),length(dL_H));
MAE_ln = zeros(length(height),length(dL_H));
MAE_beta = zeros(length(height),length(dL_H));
MRE_gh = zeros(length(height),length(dL_H));
MRE_ln = zeros(length(height),length(dL_H));
MRE_beta = zeros(length(height),length(dL_H));

%Let's also store the sigma_N_L and K_L_dB values
sigma_N_L_dB_store = zeros(length(height),length(dL_H),Num_member_UAV+1);
K_L_dB_store = zeros(length(height),length(dL_H),Num_member_UAV+1);

for j = 1:length(height)
    height(j)
    for k = 1:length(dL_H)
        dL_H(k)
        dL = sqrt(dL_H(k).^2 + height(j)^2);
        % Calculate the distances between each member UAV and the GCS
        x_i = zeros(1,Num_member_UAV); % X-coord of each member UAV, with GCS as center
        y_i = zeros(1,Num_member_UAV); % Y-coord of each member UAV, with GCS as center
        PLoS_i = zeros(1,Num_member_UAV);
        theta_i_d = zeros(1,Num_member_UAV);
        for i = 1:Num_member_UAV
            x_i(i) = dL_H(k) + radius * cos((i-1)*2*pi/Num_member_UAV);
            y_i(i) = radius * sin((i-1)*2*pi/Num_member_UAV);
            PLoS_i(i) = PLoS_v3(sqrt(x_i(i)^2+y_i(i)^2), height(j), 0, 0.1 ,7.5e-4, 8);
            theta_i_d(i) = atand(height(j)/sqrt(x_i(i)^2+y_i(i)^2));
        end
        di = sqrt(x_i.^2 + y_i.^2 + height(j)^2);
        %Calculate channel parameters affected by positions
        PLoS_L = PLoS_v3(dL_H(k),height(j),0,0.1,7.5e-4,8); % LoS Prob. btw leader & GCS
        theta_L_d = atand(height(j)/dL_H(k)); %Elevation angle of leader from GCS in degrees
        n_L = (n_min-n_max)*PLoS_L + n_max; %Path loss exponent btw leader and GCS
        n_i = (n_min-n_max).*PLoS_i + n_max; %Path loss exponent btw leader and member i (PLoS = 1)
        sigma_N_L_dB = alpha*exp(-beta*theta_L_d); sigma_N_L = 10^(sigma_N_L_dB/10);%Logarithmic Std dev of shadowing btw leader & GCS
        sigma_N_i_dB = alpha.*exp(-beta.*theta_i_d); sigma_N_i = 10.^(sigma_N_i_dB./10);%Logarithmic Std dev of shadowing btw leader & member i
        K_L_dB = K_dB_min*exp(log(K_dB_max/K_dB_min)*PLoS_L^2); K_L = 10^(K_L_dB/10); %Rician K factor btw leader and GCS
        K_i_dB = K_dB_min.*exp(log(K_dB_max/K_dB_min).*PLoS_i.^2); K_i = 10.^(K_i_dB./10); %Rician K factor btw leader and member i
        sigma_N_L_dB_store(k,j,:) = [sigma_N_L_dB, sigma_N_i_dB];
        K_L_dB_store(k,j,:) = [K_L_dB, K_i_dB];
        % Monte Carlo sim ====================================================
        v_L = sqrt(omega*K_L/(1+K_L)); s_L = sqrt(omega/(2+2*K_L)); %Rician parameters for leader UAV
        v_i = sqrt(omega.*K_i./(1+K_i)); s_i = sqrt(omega./(2+2.*K_i)); %Rician parameters for member UAV
        Pwr_sim = h_UAV*dL^(-n_L).*ricePwrrnd(v_L*ones(1, N),s_L).*10.^(normrnd(0,sigma_N_L,1,N)./10);
        I_sim = ones(length(di), N);
        for i = 1:length(di)
            I_sim(i,:) = (di(i)).^(-n_i(i)).*ricePwrrnd(v_i(i)*ones(1, N),s_i(i)).*10.^(normrnd(0,sigma_N_i(i),1,N)./10);
        end
        if length(di) > 1
            I_sim = h_UAV.*sum(I_sim);
        else
            I_sim = h_UAV.*I_sim;
        end
        SINR_sim = Pwr_sim./(I_sim+noise); %SINR
        %Get simulated SINR CDF
        mean_sinr = mean(SINR_sim); %Mean of SINR
        sigma_sinr = std(SINR_sim); %Std dev of SINR
        cdf_range = mean_sinr + 5*sigma_sinr; 
        num_bin = 150; % The number of bins to have between 0 to pdf_range
        bin_width = cdf_range / num_bin;
        h = histogram(SINR_sim,'BinWidth',bin_width,'Normalization','cdf');
        CDF_sinr_sim = h.Values(1:num_bin);
        c = conv(h.BinEdges, [0.5 0.5], 'valid'); %CDF points
        sinr_th = c(1:num_bin);
%         c_index = CDF_sinr_sim > 0;
%         sinr_th = c(c_index); 
%         CDF_sinr_sim = CDF_sinr_sim(c_index);
        %======================================================================
        % Approximation of CDF of SINR using Gauss Hermite
%         K = [K_L, K_i*ones(1,length(di))]; %Corresponding to [GCS, Member 1, Member 2]
%         sigma_N = [sigma_N_L, sigma_N_i*ones(1,length(di))]; %Corresponding to [GCS, Member 1, Member 2]
%         n = [n_L, n_i*ones(1,length(di))]; %Corresponding to [GCS, Member 1, Member 2]
        K = [K_L,K_i];
        sigma_N = [sigma_N_L,sigma_N_i];
        n = [n_L, n_i];
        hI = h_UAV * ones(1,Num_member_UAV);
        % Get the lognormal approximation CDF
        [mu_z_0,sigma_z_0] = log_normal_large_scale_approx_general(h_UAV, hI, dL, di, K(2:end), omega, sigma_N, noise_dB, n);
        [mu_z_ln,sigma_z_ln] = lognormal_rice_moment_match(eta*mu_z_0,eta*sigma_z_0,K(1),omega);
        CDF_sinr_ln = log_normal_cdf(sinr_th, eta*mu_z_ln, eta*sigma_z_ln);
        CDF_sinr_ln(isnan(CDF_sinr_ln))=0;
        err_cdf_ln = (CDF_sinr_ln-CDF_sinr_sim);
        err_cdf_ln(isnan(err_cdf_ln))=CDF_sinr_ln(isnan(err_cdf_ln));
        err_cdf_ln(isinf(err_cdf_ln))=CDF_sinr_ln(isinf(err_cdf_ln));
        MAE_ln(k,j) = mean(abs(err_cdf_ln)); %Mean abs error lognormal
        MRE_ln(k,j) = mean(abs(err_cdf_ln)./1-CDF_sinr_sim); %Mean relative error lognormal
        % Get the beta prime approximation CDF
        [k1,k2,theta1,theta2] = beta_prime_approx_v2_general(h_UAV, hI, dL, di, K, omega, sigma_N, noise_dB, n);
        CDF_sinr_beta = beta_prime_cdf(sinr_th, k1, k2, 1, theta1/theta2);
        CDF_sinr_beta(isnan(CDF_sinr_beta))=0;
        err_cdf_beta = (CDF_sinr_beta-CDF_sinr_sim);
        err_cdf_beta(isnan(err_cdf_beta))=CDF_sinr_beta(isnan(err_cdf_beta));
        err_cdf_beta(isinf(err_cdf_beta))=CDF_sinr_beta(isinf(err_cdf_beta));
        MAE_beta(k,j) = mean(abs(err_cdf_beta)); %Sum abs error beta prime
        MRE_beta(k,j) = mean(abs(err_cdf_beta)./(1-CDF_sinr_sim)); %Mean relative error beta prime
        % Get the Gauss Hermite approximation CDF
%         Np_load = load("Np_CDF_suburban_1.mat");
%         Np1 = round(interp1(Np_load.K_L_dB,Np_load.Np_converge,10*log10(K(1)))); % Num CDF terms
%         CDF_sinr_gh = gauss_hermite_rician_cdf4_v2(sinr_th, mu_z_0, sigma_z_0, K(1), omega, Np1);
%         CDF_sinr_gh(isnan(CDF_sinr_gh))=0;
%         CDF_sinr_gh(isinf(CDF_sinr_gh))=0;
%         err_cdf_gh = (CDF_sinr_gh-CDF_sinr_sim);
%         err_cdf_gh(isnan(err_cdf_gh))=CDF_sinr_gh(isnan(err_cdf_gh));
%         err_cdf_gh(isinf(err_cdf_gh))=CDF_sinr_gh(isinf(err_cdf_gh));
%         MAE_gh(k,j) = mean(abs(err_cdf_gh)); %Sum abs error Gauss Hermite
%         MRE_gh(k,j) = mean(abs(err_cdf_gh)./(1-CDF_sinr_sim)); %Mean relative error Gauss Hermite
        %======================================================================
    end
end
%%
figure()
plot(dL_H,reliability_sim);
hold on
plot(dL_H,MAE_ln,'ro');
hold on
plot(dL_H,MAE_beta,'sg');
hold on
plot(dL_H,MAE_gh,'dk');
xlabel("Distance from GCS (m)");
ylabel("Reliability");
%title("Reliability of UAV communication")
legend("Simulated","Lognormal","Beta prime","Gauss-Hermite");
%%
% Save the variable
save("MAEMRE_vs_Height_Dist_31UAV_U2G.mat",'dL_H','di','Pt','Pi','height',...
    'MAE_ln','MAE_beta','MRE_ln','MRE_beta',...
    'K_L_dB_store','sigma_N_L_dB_store');

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