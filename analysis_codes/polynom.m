function X = polynom2(d,n,z,o);

% Return polynomial expansion of data vector d to order n (each order normalised by std)
% Rik Henson 2003


if nargin<4
    o = 1;
end

if nargin<3
    z = 1;
end

d = d(:);
M = length(d);

if n>M
    error(sprintf('%dth order expansion not possible for only %d points',n,M));
end

X = []; %ones(M,1);  % 0th order (intercept)

for p=0:n
    nX = d.^p;
    if o == 1
        nX = orthog(nX,X);
    end
    if p>0 & z == 1
        nX = nX / std(nX);
    end
    X = [X nX];
end




function o=orthog(x,y);

% RH 1999

% orthog x wrt y


if(size(x,1)==1)
	x=x';
end
if(size(y,1)==1)
	y=y';
end

if isempty(y),
    o=x;    
else
    o=x-y*pinv(y)*x;
end

