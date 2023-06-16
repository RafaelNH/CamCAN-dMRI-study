function DWIsim=fun_simulate_brain_bDTI(DT,S0,data_mask,bval,bvec)

[Nx, Ny, Nz]=size(data_mask);
Nvol=length(bval);

Z6=zeros(Nvol,6);% old Z2
Ad=Z6;

DWIsim=zeros(Nx,Ny,Nz,Nvol);

for v=1:Nvol
    b=bval(v);
    Ad(v,1:6)=[b*bvec(1,v)^2, b*bvec(2,v)^2, b*bvec(3,v)^2, ...
        2*b*bvec(1,v)*bvec(2,v), 2*b*bvec(1,v)*bvec(3,v), 2*b*bvec(2,v)*bvec(3,v)];
end

%A
A=[ -Ad ones(Nvol, 1)]; %SO is now the last column

for k=1:Nz
    for j=1:Ny
        for i=1:Nx  
            if(data_mask(i,j,k)==1) 
                DTi=squeeze(DT(i,j,k,:));
                S0i=squeeze(S0(i,j,k));
                X(7)=log(S0i);
                X(1:6)=DTi;
                B=A*X';
                DWIsim(i,j,k,:)=exp(B);
            end
        end
    end
end
