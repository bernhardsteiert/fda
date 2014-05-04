% Figure 4c: IGF 20; varying AKTi
load('./Workspaces/c_signal_03302014')
extension = '03-30-2014';
addpath('./Functions/')

akti_cond = [3 24-3+1 24+3 48-3+1 48+3 7 24-7+1 24+7 48-7+1 48+7]; % IGF 20
test_cond = [1:8 17:32 41:56];
akti_dose = [.5 .1 .025 .00625 0];

colmap = jet(length(akti_dose));

figure

for itest = test_cond
    sprop = siteprop(itest,extension);
    
    if strmatch(sprop.celltype,'MCF10A','exact')
        subplot(1,2,1)
        hold on
    else
        subplot(1,2,2)
        hold on
    end
    mycol = [.7 .7 .7];
    linw = 1;
    plot(timestamp(range_ind),nanmean(c_signal(range_ind,celltype == itest),2),'Color',mycol,'LineWidth',linw)
end

legh = [];
legstr = {};
for iakti = akti_cond
    sprop = siteprop(iakti,extension);
    
    if strmatch(sprop.celltype,'MCF10A','exact')
        subplot(1,2,1)
        hold on
        title(sprintf('%s %i - MCF10A',sprop.lig_name,sprop.lig_dose))
        set(gca,'XLim',[50 400])
        set(gca,'YLim',[-.02 .03])
        ylabel('log_{10} FOXO3a [Cyt/Nuc]');
    else
        subplot(1,2,2)
        hold on
        title(sprintf('%s %i - 184A1',sprop.lig_name,sprop.lig_dose))
        set(gca,'XLim',[50 400])
        set(gca,'YLim',[-.02 .03])
        xlabel('time [min]')
    end
    mycol = colmap(akti_dose == sprop.drug_dose,:);
    linw = 2;
    legh = [legh plot(timestamp(range_ind),nanmean(c_signal(range_ind,celltype == iakti),2),'Color',mycol,'LineWidth',linw)];
    legstr{end+1} = sprintf('%s %g',sprop.drug_name,sprop.drug_dose);
    
end

legend(legh(1:length(akti_dose)),legstr(1:length(akti_dose)))