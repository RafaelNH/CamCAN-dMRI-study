% Rafael Neto Henriques, 17/03/2023
addpath('NIFTI_toolbox/')

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

nvoxel_wm=zeros(Len, 1);
nvoxel_wm_mat=zeros(Len, 48);

for s=1:Len
    % Load ROIs
    V=load_untouch_nii([pth, subject(s).name, '/ROIs/WM_JHU_ROIs.nii']);
    WM=V.img;
    
    % Load diffusion metrics
    V=load_untouch_nii([pth, subject(s).name, '/dti/FA.nii']);
    FA=V.img;
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
    MKdf_wm(s) = mean(MKdf(WM > th));
    MDdf_wm(s) = mean(MDdf(WM > th));
    OD_wm(s) = mean(OD(WM > th));
    ND_wm(s) = mean(ND(WM > th));
    Fiso_wm(s) = mean(Fiso(WM > th));
    
    
    for r=1:48
        nvoxel_wm_mat(s,r) = sum(WM(:)==r);
        if nvoxel_wm_mat(s,r)>0
            
            nvoxel_wm_mat(s, r) = sum(WM(:)== r);
            FA_wm_mat(s, r) = mean(FA(WM == r));
            MKdf_wm_mat(s, r) = mean(MKdf(WM == r));
            MDdf_wm_mat(s, r) = mean(MDdf(WM == r));
            OD_wm_mat(s, r) = mean(OD(WM == r));
            ND_wm_mat(s, r) = mean(ND(WM == r));
            Fiso_wm_mat(s, r) = mean(Fiso(WM == r));

        end
    end
    disp(s)
end

save('Diffusion_vs_Age_WMlabels_mean_noddi90', 'subject', 'sex', 'age',...
    'FA_wm','MKdf_wm','MDdf_wm','OD_wm','ND_wm','Fiso_wm',...
    'FA_wm_mat','MKdf_wm_mat','MDdf_wm_mat','OD_wm_mat','ND_wm_mat','Fiso_wm_mat',...
    'nvoxel_wm','nvoxel_wm_mat')