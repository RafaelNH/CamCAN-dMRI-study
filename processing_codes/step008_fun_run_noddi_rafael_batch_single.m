function fun_run_noddi_rafael_batch_single(pth, DWI, bvals, bvecs, mask)

addpath NODDI_toolbox_v0.9/fitting
addpath NODDI_toolbox_v0.9/models
addpath NODDI_toolbox_v0.9/models/watson
addpath niftimatlib-1.2/
addpath niftimatlib-1.2/matlab/

dwifile = [pth, DWI];
maskfile = [pth, mask];
noddi_bvals = [pth, bvals];
noddi_bvecs = [pth, bvecs];

% convert data to noddi's assumed format
noddi_data = [pth, 'NODDI_roi.mat'];
CreateROI(dwifile, maskfile, noddi_data)

% convert gradient and b-value information in noddi's assumed format
noddi_protocol = FSL2Protocol(noddi_bvals, noddi_bvecs);

% define noddi model
noddi = MakeModel('WatsonSHStickTortIsoV_B0');

% process data
noddi_params = [pth, 'NODDI_params.mat'];
batch_fitting_single(noddi_data, noddi_protocol, noddi, noddi_params);

% convert data
SaveParamsAsNIfTI(noddi_params, noddi_data, maskfile, [pth, 'NODDI'])