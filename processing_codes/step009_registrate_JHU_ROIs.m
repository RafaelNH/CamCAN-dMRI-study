function registrate_JHU_ROIs(nname,disp_mode)

reference_image=['DWI05/',nname,'/dti/FA'];
source_tc_template='/imaging/local/software/fsl/fsl64/fsl-4.1.8/fsl/data/atlases/JHU/JHU-ICBM-FA-2mm';
source_ROIs='/imaging/local/software/fsl/fsl64/fsl-4.1.8/fsl/data/atlases/JHU/JHU-WhiteMatter-labels-2mm';
output_template_dir=['DWI05/', nname,'/ROIs'];
mkdir(output_template_dir)
output_template=[output_template_dir,'/FAtemplate_in_subject_space'];
output_mat=[output_template_dir,'/flirt_par.mat'];
output_mat2=[output_template_dir,'/fnirt_par.mat'];
config_path='CamCAN/JHUFA2subject';
output_rois=[output_template_dir,'/WM_JHU_ROIs'];

if disp_mode==true
    disp('starting flirt')
end
cmd=['flirt -ref ',reference_image, ' -in ', source_tc_template,...
    ' -out ', output_template, ' -omat ', output_mat];
system(cmd)

if disp_mode==true
    disp('starting fnirt')
end

cmd=['fnirt --ref=',reference_image, ' --in=', source_tc_template,...
    ' --aff=',output_mat, ' --cout=', output_mat2, ' --config=', config_path];
system(cmd)


if disp_mode==true
    disp('starting wraping')
end

cmd=['applywarp --ref=',reference_image,' --in=', source_tc_template,...
    ' --warp=',output_mat2, ' --out=',output_template];
system(cmd)

cmd=['applywarp --ref=',reference_image,' --in=', source_ROIs,...
    ' --warp=',output_mat2, ' --out=',output_rois,' --interp=nn'];
system(cmd)

cmd=['gunzip ',output_rois,'.nii.gz'];
system(cmd)