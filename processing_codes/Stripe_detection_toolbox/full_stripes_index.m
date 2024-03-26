function fSI=full_stripes_index(SI,inputbval)

bval=load(inputbval);
MaxSI=max(SI,[],2);
SI2=SI((MaxSI<1.25),:);
Mtrend=mean(SI2,1);
fSI=SI;
ubval=unique(bval);
for b=1:length(ubval)
    SIb1=SI(:,bval==ubval(b));
    tb1=Mtrend(1,bval==ubval(b));
    TB1=[tb1; ones(1,length(tb1))]';
    Len=size(SI,1);
    for i=1:Len
        SIb1i=SIb1(i,:)';
        alfa=pinv(TB1)*SIb1i;
        E1=SIb1i-TB1*alfa;
        fSI(i,bval==ubval(b))=E1;
    end
end