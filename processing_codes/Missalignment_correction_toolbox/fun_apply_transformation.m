function fun_apply_transformation(pth, data_in, Nvol, data_out, aaa)

archstr = computer('arch');
st_comp=archstr(1:3);
if strcmp(st_comp,'win')
    sc='\';
else
    sc='/';
end

cm1=['fslsplit ', pth, sc, data_in,' ',...
    pth, sc, 'Motion_files', sc, 'iDWI',num2str(aaa)];
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
    cm2=['flirt -in ', pth, sc, 'Motion_files', sc, 'iDWI',num2str(aaa), zzz, ni,'.nii.gz',...
    ' -ref ', pth, sc, 'Motion_files', sc, 'merfin', zzz, ni, '0000.nii.gz '...
    ' -out ', pth, sc, 'Motion_files', sc, data_out, '_',zzz, ni,...
    ' -applyxfm -init ',pth, sc, 'Motion_files', sc, 'merfin',zzz, ni, '_mcf.mat/MAT_0001'];
    system(cm2);
    disp(i);
end

% save corrected DWI
cmd2=['fslmerge -t ', pth, sc, data_out];
for i=0:Nvol-1
    if i<10
        zzz='000';
    elseif i < 100
        zzz='00';
    else 
        zzz='0';
    end
    ni=num2str(i);
    cmd2=[cmd2,' ', pth, sc, 'Motion_files', sc, data_out, '_', zzz, ni,'.nii.gz'];
    disp(i);
end
system(cmd2);
cmdgz=['gunzip -f ', pth, sc, data_out, '.nii.gz'];
system(cmdgz);

cmd0=['rm -r -f ', pth, sc, 'Motion_files', sc, 'iDWI',num2str(aaa),'*'];
system(cmd0)

