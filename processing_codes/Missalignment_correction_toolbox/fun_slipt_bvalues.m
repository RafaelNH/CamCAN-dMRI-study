function fun_slipt_bvalues(selind, ...
                           sDWIname, sbvals_name, sbvecs_name, ...
                           tDWIname, tbvals_name, tbvecs_name)

% system parameters
archstr = computer('arch');
st_comp=archstr(1:3);
if strcmp(st_comp,'win')
    sc='\';
else
    sc='/';
end
addpath(['..', sc, 'NIFTI_toolbox'])

% 2 selet DWI volumes
Vs=load_untouch_nii(sDWIname);
Vs.hdr.dime.dim(5) = length(selind);
DWImat = Vs.img;
DWImat = squeeze(DWImat(:, :, :, selind));
Vs.img = DWImat;
save_untouch_nii(Vs, tDWIname);

% 3 select 
bl = load(sbvals_name);
bv = load(sbvecs_name);
bl = bl(selind);
bv = bv(:, selind);

% 4
fidl=fopen(tbvals_name,'w');
fprintf(fidl,'%4d ',bl);
fprintf(fidl,'\n');
fclose(fidl);

fidv=fopen(tbvecs_name,'w');
fprintf(fidv,'%.14f ',bv(1,:));
fprintf(fidv,'\n');
fprintf(fidv,'%.14f ',bv(2,:));
fprintf(fidv,'\n');
fprintf(fidv,'%.14f ',bv(3,:));
fprintf(fidv,'\n');
fclose(fidv);
