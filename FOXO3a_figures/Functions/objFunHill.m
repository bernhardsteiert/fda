function rssHill = objFunHill(pFit,xin,yin,qFit,pFix)
% p(1) = init; p(2) = steady-state; p(3) = IC50; p(4) = HillSlope
p = nan(1,4);
p(qFit) = pFit;
p(~qFit) = pFix;

rssHill = p(1) + (p(2)-p(1)) ./ (1 + 10.^((p(3)-xin)*p(4))) - yin;
% rssHill(1) = rssHill(1)*2;