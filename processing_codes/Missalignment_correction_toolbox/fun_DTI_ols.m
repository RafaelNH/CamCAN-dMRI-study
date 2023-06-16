function [DT,S0]=fun_DTI_ols(data_in,bval,bvec,data_mask)

%% ULLS DKI
% Implemented by Rafael Henriques
% november 2011, MRC-CBU
% adapted March 2014, MRC-CBU
%%

[Nx, Ny, Nz, Nvol]=size(data_in);


% Minimize ||AX-B||^2
% Where A, B are:

Z6=zeros(Nvol,6);% old Z2
Ad=Z6;

for v=1:Nvol
    b=bval(v);
    Ad(v,1:6)=[b*bvec(1,v)^2, b*bvec(2,v)^2, b*bvec(3,v)^2, ...
        2*b*bvec(1,v)*bvec(2,v), 2*b*bvec(1,v)*bvec(3,v), 2*b*bvec(2,v)*bvec(3,v)];
    
end

%A
A=[ -Ad ones(Nvol, 1)]; %SO is now the last column

S0=zeros(Nx,Ny,Nz);
DT=zeros(Nx,Ny,Nz,6);

for k=1:Nz
    for j=1:Ny
        for i=1:Nx  
            if(data_mask(i,j,k)==1)
                %B
                B=log(squeeze(data_in(i,j,k,:)));
                
                indinf=isinf(B);
                if sum(indinf)~=0
                    minnotinf=min(B(~indinf));
                    %maximum diffusion posible take it other
                    %direction but perhaps other aproaches
                    %can be better
                    B(indinf)=minnotinf;
                end
                
                % ULLS
                piA=pinv(A); %piA pseudoinverse of A
                X=piA*B;
                
                % diffusion tensor
                DT(i,j,k,:)=X(1:6);
                %S0
                S0(i,j,k)=exp(X(7));
                
            end
        end
    end
end
