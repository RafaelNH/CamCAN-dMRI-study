function p=fun_DWIbrain2sim_orient_general(pth,nname_s,nname_t,bvec,Nvol,step)

bv=load(bvec);
bvr=bv;

archstr = computer('arch');
st_comp=archstr(1:3);
if strcmp(st_comp,'win')
    sc='\';
else
    sc='/';
end

% Split source and target 1 and 2
cm1=['fslsplit ', nname_s,' ', pth, sc, 'Motion_files', sc, 'sou'];
system(cm1);
cm2=['fslsplit ', nname_t,' ', pth, sc, 'Motion_files', sc, 'tar'];
system(cm2);

p=zeros(Nvol-1,6);
for i=0:Nvol-1
    if i<10
        zzz='000';
    elseif i < 100
        zzz='00';
    else
        zzz='0';
    end
    ni=num2str(i);
    
    % 3 Unite each individual source and targets to use mcflirt
    cm1=['fslmerge -t ', pth, sc, 'Motion_files', sc ,'merfin', zzz, ni, ' ',...
        pth, sc, 'Motion_files', sc, 'tar', zzz, ni, '.nii.gz ',...
        pth, sc, 'Motion_files', sc, 'sou', zzz, ni, '.nii.gz'];
 
    % Unit flirt
    cm2=['mcflirt -in ', pth, sc, 'Motion_files', sc, 'merfin', zzz, ni,...
        '.nii.gz -cost normcorr -dof 6 -refvol 0 -plots -mats'];
    system(cm1);
    system(cm2);
    
    % delete slipt and merged volumes
    delete([pth, sc, 'Motion_files', sc,'merfin', zzz, ni, '.nii.gz']); % 5
    delete([pth, sc, 'Motion_files', sc,'tar', zzz, ni, '.nii.gz']); % 6
    delete([pth, sc, 'Motion_files', sc,'sou', zzz, ni, '.nii.gz']); % 7 % delete merged corrected (I cannot do this because I need the second
    % volume to be the input of the template generation iteration)
    % delete(['/imaging/rh04/camcan_dki/DWI01/', nname,'/mer', zzz, ni, '_mcf.nii.gz']);
    % Instead I will slipt it, delete the first and keep the second and latter merge it for other template generation.  
    cm1=['fslsplit ', pth, sc, 'Motion_files', sc, 'merfin', zzz, ni, '_mcf.nii.gz ',...
        pth, sc, 'Motion_files', sc,'merfin', zzz, ni];
    system(cm1);
    % now I can deleted the merged version
    delete([pth, sc, 'Motion_files', sc,'merfin', zzz, ni, '_mcf.nii.gz']);
    % delete the first
    if step==1
        delete([pth, sc, 'Motion_files', sc, 'merfin', zzz, ni, '0000.nii.gz']);
    end
    
    % motion parameters to save in one mat file 
    p0=load([pth, sc, 'Motion_files', sc, 'merfin', zzz, ni, '_mcf.par']);
    p(i+1,:)=p0(2,:);
    delete([pth, sc, 'Motion_files', sc,'merfin', zzz, ni, '_mcf.par']);
    
    % direct the bvector
    load([pth, sc, 'Motion_files', sc,'merfin', zzz, ni, '_mcf.mat',sc, 'MAT_0001'])
    bi=bv(:,i+1);
    bvr(:,i+1)=MAT_0001(1:3,1:3)*bi;
    if step==1
        rmdir([pth, sc, 'Motion_files', sc,'merfin', zzz, ni, '_mcf.mat', sc],'s')
    end
end
 
% save motion par in mat file
save([pth, sc,'motion_par_rnh'],'p');

% save new bvecs
foldout=[pth, sc,'bvecs_orient.bvec'];
fidv=fopen(foldout,'w');
fprintf(fidv,'%.14f ',bvr(1,:));
fprintf(fidv,'\n');
fprintf(fidv,'%.14f ',bvr(2,:));
fprintf(fidv,'\n');
fprintf(fidv,'%.14f ',bvr(3,:));
fprintf(fidv,'\n');
fclose(fidv);
 
% save corrected DWI
cmd2=['fslmerge -t ', pth, sc, 'DWI_aligned'];
for i=0:Nvol-1
    if i<10
        zzz='000';
    elseif i<100
        zzz='00';
    else
        zzz='0';
    end
    ni=num2str(i);
    cmd2 = [cmd2,' ', pth, sc, 'Motion_files', sc, 'merfin', zzz, ni, '0001.nii.gz'];
end
system(cmd2);

for i=0:Nvol-1
    if i<10
        zzz='000';
    elseif i<100
        zzz='00';
    else
        zzz='0';
    end
    ni=num2str(i);
    delete([pth, sc, 'Motion_files', sc,'merfin', zzz, ni, '0001.nii.gz']);
end
cmdgz=['gunzip -f ', pth, sc,'DWI_aligned.nii.gz'];
system(cmdgz);


