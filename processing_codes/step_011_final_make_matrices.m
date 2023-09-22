% Rafael Neto Henriques, 17/03/2023
addpath('/home/rh04/DKIu_v1.1/NIFTI_toolbox/')

load DWIlist.mat

Len = length(subject);

pth = 'DWI05/';
th = 0;
NODDIth = 1;
fwth = 0.9;

for s=Len:-1:1
    age(s) = subject(s).age;
    sex(s) = subject(s).sex;
    npSI(s) = subject(s).DataControl.stripes.npSIc;
end

% remove subject 548
age = age([1:547, 549:Len]);
sex = sex([1:547, 549:Len]);
npSI = npSI([1:547, 549:Len]);
subject = subject([1:547, 549:Len]);

% remove subject subjects with more than 5 volume corrupted by stripes
age = age(npSI<5);
sex = sex(npSI<5);
subject = subject(npSI<5);

Len = length(subject);

% Initialzie matrix
FA_wm_mat=zeros(Len, 48);
FA_wm=zeros(Len, 1);
MD_wm_mat=zeros(Len, 48);
MD_wm=zeros(Len, 1);
AD_wm_mat=zeros(Len, 48);
AD_wm=zeros(Len, 1);
RD_wm_mat=zeros(Len, 48);
RD_wm=zeros(Len, 1);
MKdf_wm_mat=zeros(Len, 48);
MKdf_wm=zeros(Len, 1);
MDdf_wm_mat=zeros(Len, 48);
MDdf_wm=zeros(Len, 1);
OD_wm_mat=zeros(Len, 48);
OD_wm=zeros(Len, 1);
ND_wm_mat=zeros(Len, 48);
ND_wm=zeros(Len, 1);
Fiso_wm_mat=zeros(Len, 48);
Fiso_wm = zeros(Len, 1);
MK_wm_mat=zeros(Len, 48);
MK_wm=zeros(Len, 1);
AK_wm_mat=zeros(Len, 48);
AK_wm=zeros(Len, 1);
RK_wm_mat=zeros(Len, 48);
RK_wm=zeros(Len, 1);

nvoxel_wm=zeros(Len, 1);
nvoxel_wm_mat=zeros(Len, 48);

for s=1:Len
    % Load ROIs
    V=load_untouch_nii([pth, subject(s).name, '/ROIs/WM_JHU_ROIs.nii']);
    WM=V.img;
    
    % Load diffusion metrics
    V=load_untouch_nii([pth, subject(s).name, '/dti/FA.nii']);
    FA=V.img;
    V=load_untouch_nii([pth, subject(s).name, '/dti/MD.nii']);
    MD=V.img;
    V=load_untouch_nii([pth, subject(s).name, '/dti/AD.nii']);
    AD=V.img;
    V=load_untouch_nii([pth, subject(s).name, '/dti/RD.nii']);
    RD=V.img;
    V=load_untouch_nii([pth, subject(s).name, '/dkidf/MD.nii']);
    MDdf=V.img;
    V=load_untouch_nii([pth, subject(s).name, '/dkidf/MK.nii']);
    MKdf=V.img;
    MKdf(MKdf<0)=0;
    MKdf(MKdf>5)=5;
    V=load_untouch_nii([pth, subject(s).name, '/NODDI_odi.nii']);
    OD=V.img;
    V=load_untouch_nii([pth, subject(s).name, '/NODDI_ficvf.nii']);
    ND=V.img;
    V=load_untouch_nii([pth, subject(s).name, '/NODDI_fiso.nii']);
    Fiso=V.img;
    V=load_untouch_nii([pth, subject(s).name, '/dki/MK.nii']);
    MK=V.img;
    V=load_untouch_nii([pth, subject(s).name, '/dki/AK.nii']);
    AK=V.img;
    V=load_untouch_nii([pth, subject(s).name, '/dki/RK.nii']);
    RK=V.img;
    
    if NODDIth
        WM(Fiso>fwth)=0;
    else
    % fwDTI based correction
    V=load_untouch_nii([pth, subject(s).name, '/fwdti/F.nii']);
    F=V.img;
    WM(F>fwth)=0;
    end
    
    nvoxel_wm(s) = sum(WM(:)>th);
    FA_wm(s) = mean(FA(WM > th));
    MD_wm(s) = mean(MD(WM > th));
    AD_wm(s) = mean(AD(WM > th));
    RD_wm(s) = mean(RD(WM > th));
    MKdf_wm(s) = mean(MKdf(WM > th));
    MDdf_wm(s) = mean(MDdf(WM > th));
    OD_wm(s) = mean(OD(WM > th));
    ND_wm(s) = mean(ND(WM > th));
    Fiso_wm(s) = mean(Fiso(WM > th));
    MK_wm(s) = median(MK(WM > th));
    AK_wm(s) = median(AK(WM > th));
    RK_wm(s) = median(RK(WM > th));
    
    
    for r=1:48
        nvoxel_wm_mat(s,r) = sum(WM(:)==r);
        if nvoxel_wm_mat(s,r)>0
            
            nvoxel_wm_mat(s, r) = sum(WM(:)== r);
            FA_wm_mat(s, r) = mean(FA(WM == r));
            MD_wm_mat(s, r) = mean(MD(WM == r));
            AD_wm_mat(s, r) = mean(AD(WM == r));
            RD_wm_mat(s, r) = mean(RD(WM == r));
            MKdf_wm_mat(s, r) = mean(MKdf(WM == r));
            MDdf_wm_mat(s, r) = mean(MDdf(WM == r));
            OD_wm_mat(s, r) = mean(OD(WM == r));
            ND_wm_mat(s, r) = mean(ND(WM == r));
            Fiso_wm_mat(s, r) = mean(Fiso(WM == r));
            MK_wm_mat(s, r) = median(MK(WM == r));
            AK_wm_mat(s, r) = median(AK(WM == r));
            RK_wm_mat(s, r) = median(RK(WM == r));

        end
    end
    disp(s)
end

save('Diffusion_vs_Age_WMlabels_mean_noddi90_expanded', 'subject', 'sex', 'age',...
    'FA_wm','MD_wm','AD_wm','RD_wm','MKdf_wm','MDdf_wm','OD_wm','ND_wm','Fiso_wm',...
    'MK_wm','AK_wm','RK_wm',...
    'FA_wm_mat','MD_wm_mat','AD_wm_mat','RD_wm_mat',...
    'MKdf_wm_mat','MDdf_wm_mat','OD_wm_mat','ND_wm_mat','Fiso_wm_mat',...
    'MK_wm_mat','AK_wm_mat','RK_wm_mat',...
    'nvoxel_wm','nvoxel_wm_mat')
