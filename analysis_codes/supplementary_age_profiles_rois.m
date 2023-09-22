close all
load ../data/Diffusion_vs_Age_WMlabels_mean_NODDI90_expanded.mat
load in.mat

ord_sel = 2;
pval_level = 0.05;
age = age';

% This are the representative areas that I dediced to include in
% supplementary Fig. C.1 and supplementary Fig. C.2
%roi = [17, 18]; %ALIC (Area with large age effect that maturates later)
%roi = [33, 34]; % External C (Area with large age effect that maturates later)
%roi = [45, 46]; % Uncinate fasciculus (Area that maturates later but shows lower effect with age)
%roi = [23, 24]; %Anterior Corona Radiata (Large age effect and stop maturates earlier. However, FA initial declines are likely to be confounded by ODI increases)
%roi = [43, 44]; %Sup_Fronto-Occip_Fasc (Large age effect and maturates earlier)
%roi = 3; % CCg
%roi = 4; % CCb
%roi = 5; % CCs

% Observation on other areas of interesr.
roi = [19, 20]; %PLIC (Area with early declines in FA, likely a consequence of ODI increases)
%roi = [7, 8]; %CST (Example of Area that is relatively stable with age, in line with small R2 in Fig. 3
%roi = [29, 30]; % Post_Thalamic_Radiation (region with earlier declines detected by FA, MSK, NDI)
%roi = [13, 14]; % Sup_Cerebellar_Ped - Area that shows large effects in FA/ODI - FA increases due to ODI decreases
%roi = 6; % fornix .This area is highly contaminated by free water
%increases with age. As a consequence Fiso shows large increases with age.
%As large age this region of interest basically captures free water
%(Mean(fiso)>0.8. NDI estimates for voxels containing mostly free water are
%degenerate (typically tending to the maximum value of 1) - the increase of
%NDI to 1 for this region of interest is likely to be a consequence of
%increase degeneracy of NDI estimation with age.

if length(roi) == 2
    FA_wm = (FA_wm_mat(in, roi(1)) + FA_wm_mat(in, roi(2)))/2 ;
    MDdf_wm = (MDdf_wm_mat(in, roi(1)) + MDdf_wm_mat(in, roi(2)))/2;
    MKdf_wm = (MKdf_wm_mat(in, roi(1)) + MKdf_wm_mat(in, roi(2)))/2;
    OD_wm = (OD_wm_mat(in, roi(1)) + OD_wm_mat(in, roi(2)))/2;
    ND_wm = (ND_wm_mat(in, roi(1)) + ND_wm_mat(in, roi(2)))/2;
    Fiso_wm = (Fiso_wm_mat(in, roi(1)) + Fiso_wm_mat(in, roi(2)))/2;
    nvoxel_wm = (nvoxel_wm_mat(in, roi(1)) + nvoxel_wm_mat(in, roi(2)))/2;
    AD_wm = (AD_wm_mat(in, roi(1)) + AD_wm_mat(in, roi(2)))/2;
    RD_wm = (RD_wm_mat(in, roi(1)) + RD_wm_mat(in, roi(2)))/2;
    AK_wm = (AK_wm_mat(in, roi(1)) + AK_wm_mat(in, roi(2)))/2;
    RK_wm = (RK_wm_mat(in, roi(1)) + RK_wm_mat(in, roi(2)))/2;
else
    
    FA_wm = FA_wm_mat(in, roi);
    MDdf_wm = MDdf_wm_mat(in, roi);
    MKdf_wm = MKdf_wm_mat(in, roi);
    OD_wm = OD_wm_mat(in, roi);
    ND_wm = ND_wm_mat(in, roi);
    Fiso_wm = Fiso_wm_mat(in, roi);
    nvoxel_wm = nvoxel_wm_mat(in, roi);
    AD_wm = AD_wm_mat(in, roi);
    RD_wm = RD_wm_mat(in, roi);
    AK_wm = AK_wm_mat(in, roi);
    RK_wm = RK_wm_mat(in, roi);
    
end
age = age(in);

Len = length(subject);

X = [age, age.^2];

selg1 = age >= 28 & age < 48;
selg2 = age >= 48 & age < 68;
selg3 = age >= 68 & age < 88;
xsel1 = [28, 47];
xsel2 = [48, 67];
xsel3 = [68, 87];

fig=figure;
set(fig, 'color', [1 1 1])
scatter(age, nvoxel_wm, 5, [0.4 0.4 1], 'filled')

lab2 = {'FA','MSD','MSK','NDI','ODI','F_{iso}'};


fig=figure('position', [0.0066    0.4474    1.5272    0.3184]*1e3);
set(fig, 'color', [1 1 1])
subplot(1, 6, 1)
scatter(age, FA_wm, 5, [0.4 0.4 1], 'filled')
hold on
[b, pval, r2, Age_range, yreg, yreg_inf, yreg_sup] = fun_quadratic_regression(FA_wm, age, ord_sel, pval_level);
plot(Age_range, yreg, 'black')
plot(Age_range, yreg_inf, 'black')
plot(Age_range, yreg_sup, 'black')
xlim([10 90])
if length(b) == 3
    FAm = -0.5 * b(2)/b(3);
else
    FAm = 0;
end
xlabel('Age');
title(lab2{1})
% mdl = fitlm(age(selg1), FA_wm(selg1));
% beta = mdl.Coefficients.Estimate;
% ysel1 = xsel1 * beta(2) + beta(1);
% plot(xsel1, ysel1, 'red')
%
% mdl = fitlm(age(selg2), FA_wm(selg2));
% beta = mdl.Coefficients.Estimate;
% ysel2 = xsel2 * beta(2) + beta(1);
% plot(xsel2, ysel2, 'red')
%
% mdl = fitlm(age(selg3), FA_wm(selg3));
% beta = mdl.Coefficients.Estimate;
% ysel3 = xsel3 * beta(2) + beta(1);
% plot(xsel3, ysel3, 'red')
%

%ylim([0. 1])

subplot(1, 6, 2)
scatter(age,  MDdf_wm*1000, 5, [0.4 0.4 1], 'filled')
hold on
[b, pval, r2, Age_range, yreg, yreg_inf, yreg_sup] = fun_quadratic_regression(MDdf_wm*1000, age, ord_sel, pval_level);
plot(Age_range, yreg, 'black')
plot(Age_range, yreg_inf, 'black')
plot(Age_range, yreg_sup, 'black')
xlim([10 90])
MDm = -0.5 * b(2)/b(3);
%ylim([0.4 1.4])
xlabel('Age');
title(lab2{2})

subplot(1, 6, 3)
scatter(age,  MKdf_wm, 5, [0.4 0.4 1], 'filled')
hold on
[b, pval, r2, Age_range, yreg, yreg_inf, yreg_sup] = fun_quadratic_regression(MKdf_wm, age, ord_sel, pval_level);
plot(Age_range, yreg, 'black')
plot(Age_range, yreg_inf, 'black')
plot(Age_range, yreg_sup, 'black')
xlim([10 90])
MKm = -0.5 * b(2)/b(3);
%ylim([0.4 2])
xlabel('Age');
title(lab2{3})

subplot(1, 6, 4)
scatter(age,  ND_wm, 5, [0.4 0.4 1], 'filled')
hold on
[b, pval, r2, Age_range, yreg, yreg_inf, yreg_sup] = fun_quadratic_regression(ND_wm, age, ord_sel, pval_level);
plot(Age_range, yreg, 'black')
plot(Age_range, yreg_inf, 'black')
plot(Age_range, yreg_sup, 'black')
xlim([10 90])
if length(b) == 3
    NDm = -0.5 * b(2)/b(3);
else
    NDm = 0;
end
%ylim([0.4 0.9])
xlabel('Age');
title(lab2{4})

subplot(1, 6, 5)
scatter(age,  OD_wm, 5, [0.4 0.4 1], 'filled')
hold on
[b, pval, r2, Age_range, yreg, yreg_inf, yreg_sup] = fun_quadratic_regression(OD_wm, age, ord_sel, pval_level);
plot(Age_range, yreg, 'black')
plot(Age_range, yreg_inf, 'black')
plot(Age_range, yreg_sup, 'black')
xlim([10 90])
if length(b) == 3
    ODm = -0.5 * b(2)/b(3);
else
    ODm = 0;
end
%ylim([0.0 .5])
xlabel('Age');
title(lab2{5})

subplot(1, 6, 6)
scatter(age,  Fiso_wm, 5, [0.4 0.4 1], 'filled')
hold on
[b, pval, r2, Age_range, yreg, yreg_inf, yreg_sup] = fun_quadratic_regression(Fiso_wm, age, ord_sel, pval_level);
plot(Age_range, yreg, 'black')
plot(Age_range, yreg_inf, 'black')
plot(Age_range, yreg_sup, 'black')
xlim([10 90])
%ylim([0 0.4])
if length(b) == 3
    Fim = -0.5 * b(2)/b(3);
else
    Fim = 0;
end
xlabel('Age');
title(lab2{6})
disp([FAm, MDm, MKm, NDm, ODm, Fim])

%print -f2 -depsc 'CCs.eps'

fig=figure('position', [0.0066    0.4474    1.5272    0.3184]*1e3);
set(fig, 'color', [1 1 1])
subplot(1, 4, 1)
scatter(age, AD_wm, 5, [0.4 0.4 1], 'filled')
hold on
[b, pval, r2, Age_range, yreg, yreg_inf, yreg_sup] = fun_quadratic_regression(AD_wm, age, ord_sel, pval_level);
plot(Age_range, yreg, 'black')
plot(Age_range, yreg_inf, 'black')
plot(Age_range, yreg_sup, 'black')
xlim([10 90])
if length(b) == 3
    ADm = -0.5 * b(2)/b(3);
else
    ADm = 0;
end
xlabel('Age');
title('AD')

set(fig, 'color', [1 1 1])
subplot(1, 4, 2)
scatter(age, RD_wm, 5, [0.4 0.4 1], 'filled')
hold on
[b, pval, r2, Age_range, yreg, yreg_inf, yreg_sup] = fun_quadratic_regression(RD_wm, age, ord_sel, pval_level);
plot(Age_range, yreg, 'black')
plot(Age_range, yreg_inf, 'black')
plot(Age_range, yreg_sup, 'black')
xlim([10 90])
if length(b) == 3
    RDm = -0.5 * b(2)/b(3);
else
    RDm = 0;
end
xlabel('Age');
title('RD')

subplot(1, 4, 3)
scatter(age, AK_wm, 5, [0.4 0.4 1], 'filled')
hold on
[b, pval, r2, Age_range, yreg, yreg_inf, yreg_sup] = fun_quadratic_regression(AK_wm, age, ord_sel, pval_level);
plot(Age_range, yreg, 'black')
plot(Age_range, yreg_inf, 'black')
plot(Age_range, yreg_sup, 'black')
xlim([10 90])
if length(b) == 3
    AKm = -0.5 * b(2)/b(3);
else
    AKm = 0;
end
xlabel('Age');
title('AK')

set(fig, 'color', [1 1 1])
subplot(1, 4, 4)
scatter(age, RK_wm, 5, [0.4 0.4 1], 'filled')
hold on
[b, pval, r2, Age_range, yreg, yreg_inf, yreg_sup] = fun_quadratic_regression(RK_wm, age, ord_sel, pval_level);
plot(Age_range, yreg, 'black')
plot(Age_range, yreg_inf, 'black')
plot(Age_range, yreg_sup, 'black')
xlim([10 90])
if length(b) == 3
    RKm = -0.5 * b(2)/b(3);
else
    RKm = 0;
end
xlabel('Age');
title('RK')