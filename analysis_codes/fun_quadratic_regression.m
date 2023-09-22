function [b, pval, r2, Age_range, yreg, yreg_inf, yreg_sup] = fun_quadratic_regression(y, Age, ord_sel, pval_level)

X = [Age, Age.^2, Age.^3];
Age_range = (min(Age):0.01:max(Age))';
X_range = [Age_range, Age_range.^2, Age_range.^3];

ord = ord_sel + 1;
pval = 1;
while pval > pval_level && ord>1
    ord = ord - 1;
    [b, dev, stats] = glmfit(X(:, 1:ord), y);
    pval = stats.p(ord + 1);
    disp(['order ', num2str(ord), ' pval ', num2str(pval)])
end
[b,bint,r,rint,stats] = regress(y, [ones(length(Age), 1), X(:, 1:ord)]);
yreg = [ones(length(Age_range), 1), X_range(:, 1:ord)] * b;
disp('final bs, pval and R2')
pval = stats(3);
r2 = stats(1);
MSE = stats(4);
disp(b)
disp(pval)
disp(r2)

% Convidence interval fitted values
n = length(y);
dof = n - (ord + 1);
tstudent = tinv(0.975, dof);
X = [ones(length(Age), 1), X(:, 1:ord)];
Xr = [ones(length(Age_range), 1), X_range(:, 1:ord)];
yint = zeros(length(Age_range), 1);
for i=1:length(Age_range);
    x = Xr(i, :);
    %yint(i) = tstudent*sqrt(MSE *(1 + x * pinv(X'*X) * x'));
    yint(i) = tstudent*sqrt(MSE *(x * pinv(X'*X) * x'));
end
yreg_inf = yreg - yint;
yreg_sup = yreg + yint;