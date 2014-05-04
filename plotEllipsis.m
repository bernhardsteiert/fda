function plotEllipsis(x_scores,y_scores,mycolor,alpha)

    xsize = size(x_scores);
    if xsize(1) < xsize(2)
        x_scores = x_scores';
        y_scores = y_scores';
    end

    STD = 1;                     %# standard deviations
    conf = 2*normcdf(STD)-1;     %# covers around 95% of population (for STD = 2)
    scale = chi2inv(conf,2);     %# inverse chi-squared with dof=#dimensions

    %# substract mean
    Mu = nanmean([x_scores y_scores]);
    X0 = bsxfun(@minus, [x_scores y_scores], Mu);

    %# eigen decomposition [sorted by eigen values]
    Cov = nancov(X0) * scale;
    [V D] = eig(Cov);
    [D order] = sort(diag(D), 'descend');
    D = diag(D);
    V = V(:, order);

    t = linspace(0,2*pi,100);
    e = [cos(t) ; sin(t)];        %# unit circle
    VV = V*sqrt(D);               %# scale eigenvectors
    e = bsxfun(@plus, VV*e, Mu'); %#' project circle back to orig space

    tmpx = e(1,:);
    tmpy = e(2,:);
%     ltmp = patch(tmpx, tmpy, ones(size(tmpx)), ones(size(tmpx)));
%     set(ltmp, 'FaceColor', mycolor, 'EdgeColor', 'none', 'FaceAlpha', alpha);
    plot(tmpx,tmpy,'Color',mycolor,'LineWidth',2)
    
end