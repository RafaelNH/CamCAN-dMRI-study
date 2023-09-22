close all
clear all
clc

addpath ../data
addpath NIFTI_toolbox
load Diffusion_vs_Age_WMlabels_mean_NODDI90_expanded.mat

% Labes of metrics analysed
lab = {'FA','MSD','MD', 'RD', 'AD', ...
       'MSK', 'MK', 'RK', 'AK', ...
       'NDI','ODI','F_{iso}'};

Nmes = length(lab); % number of metrics
Nsub = size(FA_wm_mat,1); % total number of subjects

age = age';
sex = sex';

% convert diffusivities to um2/ms
MDdf_wm = MDdf_wm * 1000;
MDdf_wm_mat = MDdf_wm_mat * 1000;
MD_wm = MD_wm * 1000;
MD_wm_mat = MD_wm_mat * 1000;
AD_wm = AD_wm * 1000;
AD_wm_mat = AD_wm_mat * 1000;
RD_wm = RD_wm * 1000;
RD_wm_mat = RD_wm_mat * 1000;

remove_out = 1;
% If you want to remove outliers, set remove_out to 1
%
% Note: We decided to reject outlier to minimize the influence of extrema
% values on our analysis.
% However, our analysis and be re-run without outlier exclusion by setting
% remove_out to 0. The main findings of our study maintain regardless the
% inclusion or exclusion of outliers.

% Save images - for this you need to have the JHU templates in the folder
% data and have the NIFTI_toolbox
save_images = 0;

%% Global WM
gd = [FA_wm MDdf_wm MD_wm RD_wm AD_wm MKdf_wm MK_wm RK_wm AK_wm ND_wm OD_wm Fiso_wm ];

% Save to Global WM metrics to csv file in case one want to assess then
% without running this code
fp = fopen('Global_Metrics.csv','w');
fprintf(fp,'Participant, gFA, gMSD, gMD, gRD, gAD, gMSK, gMK, gRK, gAK, gODI, gNDI, gFiso\n');
for s=1:size(gd,1)
    fprintf(fp,'%d',s);
    for m=1:Nmes
        fprintf(fp,',%8.7f',gd(s,m));
    end
    fprintf(fp,'\n');
end
fclose(fp);

zgd = zscore(gd);
%figure,boxplot(zgd)

% outliers may be due to age, so remove age (and sex) effects first
X = [polynom(age,2) sex polynom(age,2).*sex];
for m = 1:Nmes
    [t,F,p,df,R2,cR2,B,res(:,m)] = glm(gd(:,m),X,[0 1 0 0 0 0 0]',0);
end

zd = zscore(res);
figure('color', [1 1 1]), boxplot(zd)

in = [1:Nsub];
if remove_out
    out = [];
    for f = 1:size(gd,2)
        out = [out; find(abs(zd(:,f))>5)];
    end
    out  = unique(out)';
    in   = setdiff(in,out);
    Nsub = Nsub - length(out);
    
    figure('color', [1 1 1]), hist(age(out)) % any age bias?
    figure('color', [1 1 1]), boxplot(zd)
end

lab2 = {'A)  FA','B)  MSD', 'C)  MD','D)  RD','E)  AD',...
    'F)  MSK', 'G)  MK','H)  RK', 'I)  AK',...
    'J)  NDI','K)  ODI','L)  F_{iso}'};
figure(2); clf
for m = 1:Nmes
    subplot(4,3,m)
    scatter(age(in),mean(gd(in,m),2), 5, [0.4 0.4 1], 'filled')
    hold on
    rang = max(gd(in,m)) - min(gd(in,m));
    axis([15 92 min(gd(in,m)) max(gd(in,m))+0.16*rang])
    if m<4
        set(gca,'XTick',[])
    else
        xlabel('Age');
    end
    title(lab2{m}) % ylabel(lab{m})
    fprintf('\n%s\n',lab{m})
    [t,F,p,df,R2,cR2,B,res] = glm(detrend(gd(in,m),0),X(in,:),...
        [0 1 0  0  0 0 0 ]', 0); % lin
    text(15,1.1*rang+min(gd(in,m)),...
        sprintf('L:R^2=%d%%,',round(100*cR2)),'FontSize',8,...
        'Color',[1 0.33 0])
    [t,F,p,df,R2,cR2,B,res] = glm(detrend(gd(in,m),0),X(in,:),...
        [0 0 1  0  0 0 0]', 0); % quad
    text(55,1.1*rang+min(gd(in,m)),...
        sprintf('Q:R^2=%d%% ',round(100*cR2)),'FontSize',8,...
        'Color',[1 0.33 0])
    
    [t,F,p,df,R2,cR2,B,res] = glm(detrend(gd(in,m),0),X(in,:),...
        [0 0 0  1  0 0 0 ]',0); % sex
    [t,F,p,df,R2,cR2,B,res] = glm(detrend(gd(in,m),0),X(in,:),...
        [zeros(2,5) eye(2)]',0); % sex x poly
end

print -f2 -depsc 'Fig2_Scatter_Global_Age.eps'

%% ROI analyses

roi_names = textread('JHU-labels.txt','%s');

% Prepare names
roig = {}; Nroi = 0; nvox = []; new_names = {};
for r = 1:6 % unimodal ROIs (i.e. 6 first ROIs)
    
    Nroi = Nroi+1;
    roig{Nroi} = r;
    nvox(:,Nroi) = nvoxel_wm_mat(:,r);
    name = roi_names{r};
    name(find(name=='_')) = ' ';
    new_names{Nroi} = name;
    
end
for r = 1:((size(FA_wm_mat,2)-6)/2) % bimodal ROIs
    rind = (r-1)*2 + [0 1] + 6 + 1;
    tmpv = sum(nvoxel_wm_mat(:, rind), 2);
    
    Nroi = Nroi + 1;
    roig{Nroi} = rind;
    nvox(:, Nroi) = tmpv;
    name = roi_names{rind(1)}(1:(end-2));
    name(find(name == '_')) = ' ';
    new_names{Nroi} = name;
    
end

% Average bimonal ROIs and z-score
d = []; zd = [];
for r = 1:Nroi
    d  = [d; [mean(FA_wm_mat(:, roig{r}), 2),...
        mean(MDdf_wm_mat(:, roig{r}), 2),...
        mean(MD_wm_mat(:, roig{r}), 2),...
        mean(RD_wm_mat(:, roig{r}), 2),...
        mean(AD_wm_mat(:, roig{r}), 2),...
        mean(MKdf_wm_mat(:, roig{r}), 2),...
        mean(MK_wm_mat(:, roig{r}), 2),...
        mean(RK_wm_mat(:, roig{r}), 2),...
        mean(AK_wm_mat(:, roig{r}), 2),...
        mean(ND_wm_mat(:, roig{r}), 2),...
        mean(OD_wm_mat(:, roig{r}), 2),...
        mean(Fiso_wm_mat(:, roig{r}), 2)]];
    zd = [zd; [zscore(mean(FA_wm_mat(:, roig{r}), 2)),...
        zscore(mean(MDdf_wm_mat(:, roig{r}), 2)),...
        zscore(mean(MD_wm_mat(:, roig{r}), 2)),...
        zscore(mean(RD_wm_mat(:, roig{r}), 2)),...
        zscore(mean(AD_wm_mat(:, roig{r}), 2)),...
        zscore(mean(MKdf_wm_mat(:, roig{r}), 2)),...
        zscore(mean(MK_wm_mat(:, roig{r}), 2)),...
        zscore(mean(RK_wm_mat(:, roig{r}), 2)),...
        zscore(mean(AK_wm_mat(:, roig{r}), 2)),...
        zscore(mean(ND_wm_mat(:, roig{r}), 2)),...
        zscore(mean(OD_wm_mat(:, roig{r}), 2)),...
        zscore(mean(Fiso_wm_mat(:, roig{r}), 2))]];
end

fp = fopen('ROI_Metrics.csv','w');
fprintf(fp,'Participant, ROI, FA, MSD, MD, RD, AD, MSK, MK, RK, AK, NDI, ODI, Fiso \n');
for s=1:size(gd,1)
    for r=1:Nroi
        fprintf(fp,'%d,%s',s,roi_names{r});
        for m=1:Nmes
            fprintf(fp,',%8.7f',d((s-1)*Nroi+r,m));
        end
        fprintf(fp,'\n');
    end
end
fclose(fp);

Nsub = size(FA_wm_mat,1);

in = 1:Nsub;
if remove_out
    out = [];
    for m = 1:Nmes
        out = [out; find(abs(zd(:,m))>5)];
    end
    out = unique(out)';
    [sub,roi] = ind2sub([Nsub Nroi],out);
    sub = unique(sub)
    roi = unique(roi)
    in  = setdiff(in,sub);
    
    figure, hist(age(sub)) % any age bias?
    
    tmp = reshape(zd,[Nsub Nroi Nmes]);
    tmp(sub,:,:) = nan;
    tmp = reshape(tmp,[Nsub*Nroi Nmes]);
    ind = find(~isnan(tmp(:,1)));
    Nsub = Nsub - length(sub)
    
    figure, boxplot(zd(ind,:))
end

%
if save_images == 1
    v = load_untouch_nii('../data/JHU-ICBM-labels-1mm.nii');
    ROIS_img = v.img;
    v.hdr.dime.bitpix = 32;
    v.hdr.dime.datatype = 16;
    
    v0 = load_untouch_nii('../data/JHU-ICBM-FA-1mm.nii');
    SIZ = size(ROIS_img);
    R2m = zeros(SIZ);
end

% Regional Age Effects
X = polynom(age(in),2); R2 = nan(Nroi,Nmes);
rzd = reshape(zd(ind,:),[Nsub Nroi Nmes]);
rd = reshape(d(ind,:),[Nsub Nroi Nmes]);
for m = 1:Nmes
    
    for r = 1:Nroi
        [~, ~, ~, ~, R2(r,m)] = glm(squeeze(rzd(:,r,m)), X, ...
            [zeros(2,1) eye(2)]');
        
        if save_images == 1
            for sr = 1:length(roig{r})
                R2m(ROIS_img(:)==roig{r}(sr)) = R2(r, m);
            end
        end
        
    end
    if save_images == 1
        v.img = R2m;
        save_untouch_nii(v, ['../data/R2_',  lab{m},'.nii'])
    end
end
[~, rind] = sort(mean(R2, 2), 'descend');

figure(3); clf
imagesc(R2(rind, :)), colorbar, colormap(turbo(256))
set(gca, 'Xtick', 1:Nmes, 'XTickLabels', lab,...
    'Ytick', 1:27,...
    'YTickLabels', {new_names{rind}})

print -f3 -depsc -noui 'Fig3_ROIs_R2.eps'

%%

selg1 = age(in) >= 28 & age(in) < 48;
selg2 = age(in) >= 48 & age(in) < 68;
selg3 = age(in) >= 68 & age(in) < 88;

selg = [selg1, selg2, selg3];

P = nan(Nroi, Nmes, 3); L = nan(Nroi, Nmes, 3);
age_sel = age(in);
for si = 1:3
    for m = 1:Nmes
        for r = 1:Nroi
            mdl = fitlm(age_sel(selg(:, si)), squeeze(rd(selg(:,si),r,m)));
            L(r, m, si) = mdl.Coefficients.Estimate(2);
            P(r, m, si) = mdl.Coefficients.pValue(2);
        end
    end
end

Q = fdr_bh(P(:));
Q = reshape(Q, size(P));

if save_images == 1
    for si = 1:3
        for me = 1:Nmes
            Lm_pos = zeros(SIZ);
            Lm_neg = zeros(SIZ);
            for r = 1:Nroi
                if Q(r, me, si)
                    for sr = 1:length(roig{r})
                        if L(r, me, si) > 0
                            Lm_pos(ROIS_img(:)==roig{r}(sr)) = L(r,me,si);
                        else
                            Lm_neg(ROIS_img(:)==roig{r}(sr)) = -L(r,me,si);
                        end
                    end
                end
            end
            
            if max(Lm_pos(:)) > 0
                % I will save values in 1e-3
                v0.img = Lm_pos*1000;
                save_untouch_nii(v0, ['../data/L_', lab{me},...
                    '_sg', num2str(si), '_pos.nii'])
                disp(['../data/L_', lab{me},'_sg', num2str(si),'_pos.nii'])
            end
            if max(Lm_neg(:)) > 0
                v0.img = Lm_neg*1000;
                save_untouch_nii(v0, ['../data/L_',  lab{me},...
                    '_sg', num2str(si), '_neg.nii'])
                disp(['../data/L_', lab{me},'_sg', num2str(si),'_neg.nii'])
            end
        end
    end
end

%% Correlation between metrics (global)

X = polynom(age(in),2);

cc = corrcoef(zd(ind,:));
% cc = corrcoef(zgd(in,:));
cc(find(eye(size(cc)))) = 1;

% Remove
cX = kron(ones(Nroi,1),X);
azd = zd(ind,:) - cX*pinv(cX)*zd(ind,:);
%azd = zgd(in,:) - X*pinv(X)*zgd(in,:);
ac = corrcoef(azd);

lt = find(tril(ones(Nmes),-1));
cc(lt) = ac(lt);

fn = figure(4); clf
imagesc(cc),colorbar,colormap(redblue(256)),caxis([-1 1])
set(fn,'Color', [1, 1, 1])

set(gca,'Box','off');
axis off
for x = 1:Nmes
    text(x-0.2, Nmes+0.8, lab{x},'Color','k')
    text(-0.2, x, lab{x}, 'Color', 'k')
end
line([0 Nmes+0.5],[0 Nmes+0.5],'Color','k','LineWidth',2)

print -f4 -depsc 'Fig5_Metric_CorrMat.eps'


%% How many PCs? 3!
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(zd(ind,:));
%figure,bar(LATENT) % Scree criterion (visually) clearly suggests 3 comps
nPC = length(find(LATENT*length(LATENT)/sum(LATENT) > 1)) 
nPC = 4
% Kaiser criterion (normalised eigenvalues > 1) also suggests 3
fb = figure('color', [1, 1, 1]);
bar(EXPLAINED, 'FaceColor', [0.4 0.4 1]);
ax = gca;
ax.FontSize = 20; 
set(gca,'XTick',2:2:12);
eval(sprintf('print -depsc -f%d Fig_explained_variance.eps',...
    fb.Number))
EXPLAINED
sum(EXPLAINED(1:nPC)) % >95% 

%% Factor analysis (varimax rotation)
[LAMBDA, PSI, T, STATS, F] = factoran(zd(ind,:), nPC); 
%, 'Rotate', 'promax'); Promax similar to default varimax

f6 = figure('color', [1 1 1],...
    'position', [599.4, 126.6, 721.6, 635.4]); clf
for f = 1:4
    subplot(6, 1, f),
    b = bar(LAMBDA(:, f));
    b.FaceColor = [1 0.33 0];
    b.EdgeColor = [0 0 0];
    set(gca,'XTickLabel',lab)
    axis([0 Nmes+1 -1.2 1.2])
    title(sprintf('Factor %d',f))
end
% 1st PC ~= +FA, -MD, +MK, +ND
% 2nd PC ~= +MD, +Fiso
% 3rd PC ~= -FA, +OD

f2 = figure(8); %set(f2,'Position',[500 500 1500 1500]); 
clf
f3 = figure(9); %set(f3,'Position',[1000 1000 1500 1500]);

f2 = figure(8); set(f2,'Position',[100 30 850 740]);
clf
f3 = figure(9); set(f3,'Position',[100 30 850 740]);
clf
for f = 1:nPC
    tmp = zd(ind,:) * LAMBDA(:,f);
    tmp = reshape(tmp,[Nsub Nroi]);
    
    figure(f6)
    subplot(3, nPC, f+2*nPC), 
    scatter(age(in), mean(tmp,2), 5, [0.4 0.4 1], 'filled')
    title(sprintf('Factor %d',f))
    xlabel('Age'), 
    ylabel('Average Loading across ROIs')
    axis([min(age) max(age) -6 6]);
    
    % linear
    [t,F,p,df,R2,cR2,B,res] = glm(detrend(mean(tmp,2),0), X, [0 1 0]', 0);
    text(20, 5.5, sprintf('L:R^2=%d%%,',round(100*cR2)), ...
        'FontSize', 8, 'Color', [1 0.33 0])
    
    % quadratic
    [t,F,p,df,R2,cR2,B,res] = glm(detrend(mean(tmp,2),0), X, [0 0 1]', 0);
    text(55, 5.5, sprintf('Q:R^2=%d%% ',round(100*cR2)),...
        'FontSize', 8, 'Color', [1 0.33 0])
    
    figure(f2), subplot(5,1,f),hold on
    mb = mean(tmp,1); sb = std(tmp)/sqrt(Nsub);
    mb = detrend(mb,0); % absolute value has little meaning?
    b = bar(mb);
    b.FaceColor = [1 0.33 0];
    b.EdgeColor = [0 0 0];
    title(sprintf('ROI Loadings: Factor %d',f))
    axis([0.5 27.5 -0.2 0.2]);
    set(gca,'XTick',[],'YTick',[-0.1 0 0.1]);
    
    figure(f3), subplot(5,1,f),hold on
    title(sprintf('Effect of Age (R^2): Factor %d',f)) 
    for r = 1:Nroi
        [~,~,~,~,cb(r)] = glm(tmp(:,r),X,[zeros(2,1) eye(2)]');
    end
    b = bar(cb);
    b.FaceColor = [1 0.33 0];
    b.EdgeColor = [0 0 0];
    axis([0.5 25.5 0 0.7]);
    set(gca,'XTick',[],'YTick',[0 0.5]);
end

eval(sprintf('print -depsc -f%d Fig6_Factors_Scores_and_age.eps',...
    f6.Number))

figure(f2), hold on
subplot(4,1,4),axis([0.5 27.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
set(gca,'color','none','XTick',[],'YTick',[],'XColor', 'none',...
    'YColor','none');
eval(sprintf('print -depsc -f%d Fig7_ROI_loadings.eps',f2.Number))

figure(f3), hold on
subplot(4,1,4),axis([0.5 27.5 0 1])
for r = 1:Nroi
    t = text(r,0,new_names{r},'Rotation',90,'FontSize',10);
end
set(gca,'color','none','XTick',[],'YTick',[],'XColor', 'none', ...
    'YColor','none');
eval(sprintf('print -depsc -f%d Fig8_ROI_Age_Cor.eps',f3.Number))

