function area = ellipsis_area(x_scores,y_scores)

xsize = size(x_scores);
if xsize(1) < xsize(2)
    x_scores = x_scores';
    y_scores = y_scores';
end

magic_factor = 2;
y_scores = y_scores * magic_factor; 

STD = 1;                     %# standard deviations
conf = 2*normcdf(STD)-1;     %# covers around 95% of population (for STD = 2)
scale = chi2inv(conf,2);     %# inverse chi-squared with dof=#dimensions

%# substract mean
Mu = mean([x_scores y_scores]);
X0 = bsxfun(@minus, [x_scores y_scores], Mu);

%# eigen decomposition [sorted by eigen values]
Cov = cov(X0) * scale;
[V D] = eig(Cov);
D = sort(diag(D), 'descend');
Dsqrt = sqrt(D);

area = pi*Dsqrt(1)*Dsqrt(2);

end