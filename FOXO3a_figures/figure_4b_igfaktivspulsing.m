% Figure 4b: IGF + AKT compared to pulsing conditions
load('./Workspaces/c_signal_03302014')
extension = '03-30-2014';
addpath('./Functions/')

sites_all = 1:72;
extension = '03-30-2014';

puls_thres = .55;

ratio_puls = nan(size(sites_all));
ratio_puls_bd = nan(2,length(sites_all));
rng(0) % Make sure that bootstrap samples are reproducible
ratio_fun = @(x) sum(x > puls_thres) / length(x);

for isite = sites_all
    ratio_puls(isite) = sum(dists(celltype == isite) > puls_thres) / sum(celltype == isite);
    ratio_puls_bd(:,isite) = bootci(2000,{ratio_fun,dists(celltype == isite)},'alpha',.32);
end

reference_cond = [65:67 69:71];
test_cond = [1:8 17:32 41:56];

chi2s = nan(length(reference_cond),length(test_cond));
for iref = reference_cond
    sref = siteprop(iref,extension);
    if strmatch(sref.celltype,'MCF10A','exact')
        mycol = [0 0 0];
    else
        mycol = [1 0 0];
    end
    
    for itest = test_cond
        stest = siteprop(itest,extension);
        if strmatch(sref.celltype,stest.celltype,'exact')
            chi2s(iref==reference_cond,itest==test_cond) = nansum((nanmean(c_signal(range_ind,celltype == iref),2)-nanmean(c_signal(range_ind,celltype == itest),2)).^2);
        end
    end
    
end

[tmp opt_cond] = min(chi2s,[],2);

reference_cond = [65:67 69:71];
test_cond = [1:8 17:32 41:56];
unstim_cond = [68 72];

f10 = figure;

xfac = 1;
yfac = .6;
fontsize = 8;

setFigure(f10,xfac,yfac,fontsize)

mysb = [1:3 6:8];
for isb = 1:length(mysb)
    subplot(2,5,mysb(isb))
    hold on
    legh = [];
    legstr = {};
    
    sref = siteprop(reference_cond(isb),extension);
    for itest = test_cond
        stest = siteprop(itest,extension);
        if strmatch(sref.celltype,stest.celltype,'exact')
            if itest == test_cond(opt_cond(isb))
                mycol = [1 0 0];
                linw = 2;
                legh = [legh plot(timestamp(range_ind),nanmean(c_signal(range_ind,celltype == itest),2),'Color',mycol,'LineWidth',linw)];
                legstr{end+1} = sprintf('%s %g; %s %g',stest.lig_name,stest.lig_dose,stest.drug_name,stest.drug_dose);
            else
                mycol = [.7 .7 .7];
                linw = 1;
                plot(timestamp(range_ind),nanmean(c_signal(range_ind,celltype == itest),2),'Color',mycol,'LineWidth',linw)
            end
        end
    end
    legh = [legh plot(timestamp(range_ind),nanmean(c_signal(range_ind,celltype == reference_cond(isb)),2),'Color',[0 0 1],'LineWidth',2)];
    legstr{end+1} = sprintf('%s %g',sref.lig_name,sref.lig_dose);
    legend(legh,legstr,'FontSize',6)
    
    for iunstim = unstim_cond
        sunstim = siteprop(iunstim,extension);
        if strmatch(sref.celltype,sunstim.celltype,'exact')
            plot(timestamp(range_ind),nanmean(c_signal(range_ind,celltype == iunstim),2),'Color',[1 0 1],'LineWidth',2)
        end
    end
    
    set(gca,'XLim',[50 400],'XTick',70:50:370,'XTickLabel',-50:50:250)
    set(gca,'YLim',[-.02 .03])
    
    title(sprintf('%i: %s %i - %s',isb,sref.lig_name,sref.lig_dose,sref.celltype))
%     title(sprintf('Condition %i - %s',isb,sref.celltype))
    
    if mysb(isb) == 1
        ylabel('log 10 FOXO ratio')
    end
    if mysb(isb) == 6
        xlabel('time [min]')
    end
end

subplot(2,5,[4 5 9 10])
h = bar([ratio_puls(reference_cond); ratio_puls(test_cond(opt_cond)); ratio_puls(unstim_cond([1 1 1 2 2 2]))]');
set(h(1),'FaceColor',[0 0 1])
set(h(2),'FaceColor',[1 0 0])
set(h(3),'FaceColor',[1 0 1])
legend('Control','IGF / AKT combination','Unstimulated','Location','NorthWest')
ylabel('ratio pulsing cells')
set(gca,'XLim',[.5 6.5])
x1 = get(get(h(1),'children'), 'xdata');
x1 = mean(x1([1 3],:));
x2 = get(get(h(2),'children'), 'xdata');
x2 = mean(x2([1 3],:));
x3 = get(get(h(3),'children'), 'xdata');
x3 = mean(x3([1 3],:));

hold on
h=errorbar(x1,ratio_puls(reference_cond),ratio_puls(reference_cond)-ratio_puls_bd(1,reference_cond),ratio_puls_bd(2,reference_cond)-ratio_puls(reference_cond),'k'); set(h,'linestyle','none')
h=errorbar(x2,ratio_puls(test_cond(opt_cond)),ratio_puls(test_cond(opt_cond))-ratio_puls_bd(1,test_cond(opt_cond)),ratio_puls_bd(2,test_cond(opt_cond))-ratio_puls(test_cond(opt_cond)),'k'); set(h,'linestyle','none')
h=errorbar(x3,ratio_puls(unstim_cond([1 1 1 2 2 2])),ratio_puls(unstim_cond([1 1 1 2 2 2]))-ratio_puls_bd(1,unstim_cond([1 1 1 2 2 2])),ratio_puls_bd(2,unstim_cond([1 1 1 2 2 2]))-ratio_puls(unstim_cond([1 1 1 2 2 2])),'k'); set(h,'linestyle','none')
