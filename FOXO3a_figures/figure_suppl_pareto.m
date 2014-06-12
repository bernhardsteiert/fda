% Figure supplementary: %variance explained vs. #basis functions
addpath('./Functions/')

load('./Workspaces/harm_basis_fPCA_5basis_noFGF')

f1 = figure;

xfac = 1;
yfac = 1;
fontsize = 24;

setFigure(f1,xfac,yfac,fontsize)

thres_var = 0.95;

propvar = squeeze(sum(harmscr.^2));
propvar = propvar./sum(propvar);
propvar = propvar.*sum(c_signal_pcastr.varprop(1:length(propvar)));

cumprobs = cumsum([0;propvar']);
[tmp thres_ind] = min(abs(cumprobs - thres_var));
if cumprobs(thres_ind)-thres_var < 0
    thres_ind = thres_ind + 1;
end

dx = (thres_var - cumprobs(thres_ind)) ./ propvar(thres_ind-1);

plot(0:length(propvar),cumprobs,'ko-')
hold on
plot([0 thres_ind-1+dx],[1 1]*thres_var,'r--')
plot([thres_ind-1+dx thres_ind-1+dx],[0 thres_var],'r--')
set(gca,'XLim',[0 5])
% set(gca,'YScale','log')

xlabel('fPCA basis functions')
ylabel('cumulative variance explained')
