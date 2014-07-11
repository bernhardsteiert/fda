% Figure 5b: EGF vs. MEKi heatmap
addpath('./Functions/')
load('Workspaces/scores_04182014')
extension = '04-18-2014';
sites_all = [1:39 41:72];

pcs = [2 3];

figure

egf_dose = [100 20 4 .8 .16 0];
meki_dose = [.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0];
myct = {'MCF10A','184A1'};

mean_pcs = nan(length(egf_dose),length(meki_dose),length(pcs),2); % IGF / AKTi / PCs / Celltype

for isite = sites_all
    s = siteprop(isite,extension);
    
    Mu = nanmean(scores_all(pcs,celltype == isite),2);
    mean_pcs(s.lig_dose == egf_dose, s.drug_dose == meki_dose, 1, strmatch(s.celltype,myct,'exact')) = Mu(1);
    mean_pcs(s.lig_dose == egf_dose, s.drug_dose == meki_dose, 2, strmatch(s.celltype,myct,'exact')) = Mu(2);

end

for icondpc = 1:size(mean_pcs,3)
    for icondct = 1:size(mean_pcs,4)
        subplot(size(mean_pcs,3),size(mean_pcs,4),(icondpc-1)*size(mean_pcs,4) + icondct)
        imagesc(mean_pcs(:,:,icondpc,icondct))
        title(sprintf('%s - PC %i',myct{icondct},pcs(icondpc)))
        colorbar
        set(gca,'XTick',1:length(meki_dose),'XTickLabel',meki_dose)
        set(gca,'YTick',1:length(egf_dose),'YTickLabel',egf_dose)
        xlabel('MEKi [muM]')
        ylabel('EGF [ng/mL]')
        
        colormap(flipud(hot))
    end
end