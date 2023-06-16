function apply_trans_other_data(pth, data2trans, Nvol)

archstr = computer('arch');
st_comp=archstr(1:3);
if strcmp(st_comp,'win')
    sc='\';
else
    sc='/';
end

cm1=['fslsplit ', pth, sc, data2trans,' ',...
    pth, sc, 'Motion_files', sc, 'iDWI'];
system(cm1);
for i=0:Nvol-1
    if i<10
        zzz='000';
    elseif i < 100
        zzz='00';
    else
        zzz='0';
    end
    ni=num2str(i);
    cm2=['flirt -in ', pth, sc, 'Motion_files', sc, 'iDWI',zzz, ni,'.nii.gz',...
    ' -ref ', pth, sc, 'Motion_files', sc, 'merfin', zzz, ni, '0000.nii.gz '...
    ' -out ', pth, sc, 'Motion_files', sc, 'DWI_motion_',zzz, ni,...
    ' -applyxfm -init ',pth, sc, 'Motion_files', sc, 'merfin',zzz, ni, '_mcf.mat/MAT_0001'];
    system(cm2);
end

% save corrected DWI
cmd2=['fslmerge -t ', pth, sc, 'DWI_motion'];
for i=0:Nvol-1
    if i<10
        zzz='000';
    elseif i < 100
        zzz='00';
    else 
        zzz='0';
    end
    ni=num2str(i);
    cmd2=[cmd2,' ', pth, sc, 'Motion_files', sc,'DWI_motion_', zzz, ni,'.nii.gz'];
end
system(cmd2);
cmdgz=['gunzip -f ', pth, sc,'DWI_motion.nii.gz'];
system(cmdgz);

