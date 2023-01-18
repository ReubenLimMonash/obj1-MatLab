% Date: 6/1/2023
% Desc: To test the convergence of CDF at 10^8 samples, whether the values
% are consistent up to two decimal places.

cdf_sim_10e8 = load("CDF_vs_Height_Dist_10e8/CDF_vs_Height_Dist_31UAV.mat", "CDF_sim").CDF_sim;
cdf_sim_10e7 = load("CDF_vs_Height_Dist_10e7/CDF_vs_Height_Dist_31UAV_10e7.mat", "CDF_sim").CDF_sim;

diff_cdf_sim = abs(floor(100.*cdf_sim_10e8)-floor(100.*cdf_sim_10e7));
max_diff_cdf_sim = max(diff_cdf_sim,[],"all");
conv_sim = sum(diff_cdf_sim > 0,'all');
percent_conv = (numel(diff_cdf_sim)-conv_sim)/numel(diff_cdf_sim);

diff_cdf_sim_rounded = abs(floor(100.*round(cdf_sim_10e8,2))-floor(100.*round(cdf_sim_10e7,2)));
max_diff_cdf_sim_rounded = max(diff_cdf_sim_rounded,[],"all");
conv_sim_rounded = sum(diff_cdf_sim_rounded > 0,'all');
percent_conv_rounded = (numel(diff_cdf_sim_rounded)-conv_sim_rounded)/numel(diff_cdf_sim_rounded);

%%

dir_10e8 = "CDF_vs_Height_Dist_10e8";
dir_10e7 = "CDF_vs_Height_Dist_10e7";

files_10e8 = dir(fullfile("CDF_vs_Height_Dist_10e8","*.mat"));

min_percent_conv = 100;
max_maxdiff =0;

for i = 1:length(files_10e8)
    file_name = files_10e8(i).name;
    file_10e8 = dir_10e8 + "/" + file_name;
    file_10e7 = dir_10e7 + "/" + extractBefore(file_name,length(file_name)-3) + "_10e7.mat";
    if isfile(file_10e7)
        cdf_sim_10e8 = load(file_10e8, "CDF_sim").CDF_sim;
        cdf_sim_10e7 = load(file_10e7, "CDF_sim").CDF_sim;
        diff_cdf_sim = abs(floor(10.*cdf_sim_10e8)-floor(10.*cdf_sim_10e7));
        max_diff_cdf_sim = max(diff_cdf_sim,[],"all");
        conv_sim = sum(diff_cdf_sim > 0,'all');
        percent_conv = (numel(diff_cdf_sim)-conv_sim)/numel(diff_cdf_sim);
        if percent_conv < min_percent_conv
            min_percent_conv = percent_conv
            min_percent_conv_file = file_10e8
        end
        if max_diff_cdf_sim > max_maxdiff
            max_maxdiff = max_diff_cdf_sim
            max_maxdiff_file = file_10e8
        end
    else
        alert = "File " + file_10e7 + " not found!";
    end
end