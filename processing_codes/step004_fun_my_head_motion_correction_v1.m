function P1 = fun_my_head_motion_correction_v1(pth, selind, DWIname_pre, ... 
                                             bvals, bvecs, mask, ...
                                             DWIname_orig)
% pth - path where images are saved
% selind - index of volumes where motion did not occour (FSL index + 1)
% DWIname_pre - Name of the DWI file skull scripted (nii)
% bvals - Name of the bvals file
% bvecs - Name of the bvecs file
% DWIname_orig - Original data to be corrected

disp('Setting paths and pre-processing files ...')
% 0.1) define system parameters and add paths
archstr = computer('arch');
st_comp=archstr(1:3);
if strcmp(st_comp,'win')
    sc='\';
else
    sc='/';
end
addpath(['..', sc, 'NIFTI_toolbox'])

% 0.2) make output folder
mkdir([pth, sc, 'Motion_files'])

% 0.3) generate a smooth version of the data. I notice that the procedures
% is more stable if we smoothing the data before generating the templates
% and computing the motion parameters. Note that this no mean that I am
% correcting a smoothed version of the data, at the end motion parameters
% are applied on the data pointed in output DWIname_orig 
V=load_untouch_nii([pth, sc, DWIname_pre]);
DWI=double(V.img);
sd = 1.5 /(sqrt(8*log(2))); % Convert FWHM to sd
S = size(DWI);
Nvol = S(4);
size_k=[7 7 7];
for im=Nvol:-1:1
    A = squeeze(DWI(:, :, :, im));
    W = smooth3(A, 'gaussian', size_k, sd);
    DWI(:,:,:,im)=W;
end
DWI_smooth_name = [pth, sc, 'Motion_files', sc, 'DWI_GF.nii'];
V.img=DWI;
save_untouch_nii(V, DWI_smooth_name)

disp('Selection of volumes to generate the template ...')

% 1) Selection of volumes to generate the templates
bvals = [pth, sc, bvals];
bvecs = [pth, sc, bvecs];
DWIname_sel = [pth, sc, 'Motion_files', sc, 'DWI_sel.nii'];
bvals_sel = [pth, sc, 'Motion_files', sc, 'bvals_sel.bval'];
bvecs_sel = [pth, sc, 'Motion_files', sc, 'bvecs_sel.bvec'];
fun_slipt_bvalues(selind, ...
                  DWI_smooth_name, bvals, bvecs, ...
                  DWIname_sel, bvals_sel, bvecs_sel)

disp('Generate the template (preliminary version)...')

% 2) Compute DTI (I am using the simple ols DTI approach only to reduce
% time, more accurate fitting procedures should not be necessary, since
% I only want to generate a plausible contrast)
V_in = load_untouch_nii(DWIname_sel);
data_in = double(V_in.img);
bval = load(bvals_sel);
bvec = load(bvecs_sel);
mask = [pth, sc, mask];
V_mask = load_untouch_nii(mask);
data_mask = double(V_mask.img);
[DT, S0] = fun_DTI_ols(data_in, bval, bvec, data_mask);

% 3) Generate templates based on the selected data. At this point I only
% generate templates for the selected data (see reason of the next step)
DWIsim = fun_simulate_brain_bDTI(DT, S0, data_mask, bval, bvec);
Vol_sim=V_in;
Vol_sim.hdr.dime.dim(5)=length(selind);
Vol_sim.hdr.dime.pixdim(5)=0;
Vol_sim.hdr.dime.xyz_t=0;
Vol_sim.hdr.dime.datatype=16;
Vol_sim.img=DWIsim;
template1 = [pth, sc, 'Motion_files',sc,'DWI_simDTI_1st.nii'];
save_untouch_nii(Vol_sim, template1);

disp('Ajust template (final version)...')
% 4) Correct volumes seleted - this is just an attempt to correct some residual
% motion on the selected volumes)
p0 = fun_DWIbrain2sim_orient_general(pth, DWIname_sel, template1, bvecs_sel, length(selind), 1);

P1(:,:,1)=p0;

%% Step 2 - simulate templates for all the data
% Produce the final version of the diffusion tensor 
foldbvec = [pth, sc, 'bvecs_orient.bvec'];
DWI = [pth, sc,'DWI_aligned.nii'];
V_in = load_untouch_nii(DWI);
data_in = double(V_in.img);
bvec = load(foldbvec);
[DT,S0] = fun_DTI_ols(data_in, bval, bvec, data_mask);


% Update the bvectors that were realign in the previous step
foldbvec=[pth, sc,'bvecs_orient.bvec'];
foldbvec_all=[pth, sc,'bvecs.bvec'];
foldbval=[pth, sc,'bvals.bval'];
bval=load(foldbval);
b1000o=load(foldbvec);
ball=load(foldbvec_all);

ball(:, selind) = b1000o;

fidv=fopen([pth, sc,'bvecs_orient.bvec'],'w');
fprintf(fidv,'%.14f ',ball(1,:));
fprintf(fidv,'\n');
fprintf(fidv,'%.14f ',ball(2,:));
fprintf(fidv,'\n');
fprintf(fidv,'%.14f ',ball(3,:));
fprintf(fidv,'\n');
fclose(fidv);

% Produce the last version of the templates (now all volumes are taken into
% account)
DWIsim = fun_simulate_brain_bDTI(DT,S0,data_mask,bval,ball);
Vol_sim=V_in;
[Nx, Ny, Nz, Nv] = size(DWIsim);
Vol_sim.hdr.dime.dim(5)=Nv;
Vol_sim.hdr.dime.pixdim(5)=0;
Vol_sim.hdr.dime.xyz_t=0;
Vol_sim.hdr.dime.datatype=16;
Vol_sim.img=DWIsim;
save_untouch_nii(Vol_sim,[pth, sc, 'Motion_files', sc,'ALLDWIsim_brain_bDTI.nii']);

disp('Compute head motion parameters ...')

% Use templates to correct the volumes
foldbvec=[pth, sc,'bvecs_orient.bvec'];
nname_t=[pth, sc, 'Motion_files',sc, 'ALLDWIsim_brain_bDTI.nii'];
nname_s=[pth, sc, 'Motion_files',sc, 'DWI_GF.nii'];
P1=fun_DWIbrain2sim_orient_general(pth,nname_s,nname_t,foldbvec,length(bval), 2);

% Apply correction to data given in input DWIname_orig
disp('Correct Head motion missalignments ...')
apply_trans_other_data(pth, DWIname_orig, length(bval));

% delete some files that you should not need
cmd =[ 'rm -f -r ', pth, sc, 'Motion_files'];
system(cmd);

