% Figure 5b: EGF vs. MEKi fPCA
addpath('./Functions/')
load('Workspaces/scores_04182014')
extension = '04-18-2014';
sites_all = [1:39 41:72];

pcs = [2 3];

figure
hold on

fixedEGF = 13:24;
fixedMEKi = [4 24-4+1 24+4 48-4+1 48+4 72-4+1 9 24-9+1 24+9 48-9+1 48+9 72-9+1];
fixedEGF = fixedEGF([6:-1:1 7:12]);
fixedMEKi = fixedMEKi([6:-1:1 7:12]);

ligand_dose = [100 20 4 .8 .16 0];
drug_dose = [.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0];

colmap = jet(length(drug_dose));

for iplot = [1 2 4 5]
    subplot(2,3,iplot)
    plot(scores_all(pcs(1),:),scores_all(pcs(2),:),'.','Color',[.7 .7 .7])
    set(gca,'XLim',[-.25 .25],'YLim',[-.15 .1])
end

legstr = cell(1,length(drug_dose));
for isite = fixedEGF
    sprop = siteprop(isite,extension);
    if strmatch(sprop.celltype,'MCF10A','exact')
        subplot(2,3,1)
        hold on
        legstr{sprop.drug_dose == drug_dose} = sprintf('MEKi %0.2f',1000*drug_dose(sprop.drug_dose == drug_dose));
    else
        subplot(2,3,2)
        hold on
    end
    
    plot(scores_all(pcs(1),celltype == isite),scores_all(pcs(2),celltype == isite),'o','Color',colmap(sprop.drug_dose == drug_dose,:),'MarkerFaceColor',colmap(sprop.drug_dose == drug_dose,:))
    plotEllipsis(scores_all(pcs(1),celltype == isite),scores_all(pcs(2),celltype == isite),colmap(sprop.drug_dose == drug_dose,:),.5);
    title(sprintf('%s %s %g ng/ml; varying %s',sprop.celltype,sprop.lig_name,sprop.lig_dose,sprop.drug_name))
end

subplot(2,3,3)
set(gca,'Visible','off')
set(gca,'CLim',[0 1])
colormap(colmap)
colorbar('Location','West','YTick',linspace(1./(2*length(drug_dose)),1-1./(2*length(drug_dose)),length(drug_dose)),'YTickLabel',legstr,'TickLength', [0 0]) % Vertical colorbar

legstr = cell(1,length(drug_dose));
for isite = fixedMEKi
    sprop = siteprop(isite,extension);
    if strmatch(sprop.celltype,'MCF10A','exact')
        subplot(2,3,4)
        hold on
        legstr{sprop.lig_dose == ligand_dose} = sprintf('EGF %g',ligand_dose(sprop.lig_dose == ligand_dose));
    else
        subplot(2,3,5)
        hold on
    end
    
    try
        plot(scores_all(pcs(1),celltype == isite),scores_all(pcs(2),celltype == isite),'o','Color',colmap(sprop.lig_dose == ligand_dose,:),'MarkerFaceColor',colmap(sprop.lig_dose == ligand_dose,:))
        plotEllipsis(scores_all(pcs(1),celltype == isite),scores_all(pcs(2),celltype == isite),colmap(sprop.lig_dose == ligand_dose,:),.5);
        title(sprintf('%s - MEKi 1/640; varying EGF',sprop.celltype))
    end
end

subplot(2,3,6)
set(gca,'Visible','off')
set(gca,'CLim',[0 1])
colormap(colmap)
colorbar('Location','West','YTick',linspace(1./(2*length(ligand_dose)),1-1./(2*length(ligand_dose)),length(ligand_dose)),'YTickLabel',legstr,'TickLength', [0 0]) % Vertical colorbar