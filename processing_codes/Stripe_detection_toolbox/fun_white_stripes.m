function Index_white=fun_white_stripes(inputDWIfold,inputbval)

V_DWI=spm_vol(inputDWIfold); DWI_data=spm_read_vols(V_DWI);

bval=load(inputbval);

siz=size(DWI_data);
kn=[0 -1 2 -1 0; 0 -1 2 -1 0; 0 -1 2 -1 0; 0 -1 2 -1 0; 0 -1 2 -1 0];
Data_lines=zeros(siz(1),siz(2)+4,siz(3)+4,siz(4));
Index_white=zeros(1,siz(4));
for vol=1:siz(4)
    for sag=1:siz(1)
    sag_image=squeeze(DWI_data(sag,:,:,vol));
    sag_lines=conv2(sag_image,kn);
    Data_lines(sag,:,:,vol)=sag_lines;
    end
    Index_white(vol)=sum(sum(sum(abs(squeeze(Data_lines(:,:,:,vol))))));
end

bval = round(bval / 100)*100;
ubval=unique(bval);

for b=1:length(ubval)
m=min(Index_white(bval==ubval(b)));
Index_white(bval==ubval(b))=Index_white(bval==ubval(b))/m;
end
