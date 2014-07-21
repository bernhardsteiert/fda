% Figure 5c: Bar charts EGF4 varying MEKi and IGF varying AKTi
egf_dose = 4; % Choose EGF dose here
igf_dose = 20; % Same for IGF

addpath('./Functions/')
egfmeki = load('Workspaces/dists_04182014');
igfakti = load('./Workspaces/c_signal_03302014');

sites_egfmeki = [1:39 41:72];
sites_igfakti = 1:72;

extension_egfmeki = '04-18-2014';
extension_igfakti = '03-30-2014';

puls_thres = .3;

meki_doses = [.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0];
akti_doses = [.5 .1 .025 .00625 0];

egfmeki_pulsing = nan(1,length(meki_doses)+1,2);
igfakti_pulsing = nan(1,length(akti_doses)+1,2);

egfmeki_pulsing_bd = nan(2,length(meki_doses)+1,2);
igfakti_pulsing_bd = nan(2,length(akti_doses)+1,2);

rng(0) % Make sure that bootstrap samples are reproducible
ratio_fun = @(x) sum(x > puls_thres) / length(x);

for i = 1:length(sites_egfmeki)
    
    isite = sites_egfmeki(i);
    s = siteprop(isite,extension_egfmeki);
    
    icell = 0;
    if ~isempty(strmatch(s.celltype,'MCF10A','exact'))
        icell = 1;
    elseif ~isempty(strmatch(s.celltype,'184A1','exact'))
        icell = 2;
    end
    
    if s.lig_dose == egf_dose && icell > 0  
        egfmeki_pulsing(1,s.drug_dose == meki_doses,icell) = sum(egfmeki.dists(egfmeki.celltype == isite) > puls_thres) / sum(egfmeki.celltype == isite);
        egfmeki_pulsing_bd(:,s.drug_dose == meki_doses,icell) = bootci(2000,{ratio_fun,egfmeki.dists(egfmeki.celltype == isite)},'alpha',.32);
    elseif s.lig_dose == 0 && s.drug_dose == 0 && icell > 0  
        egfmeki_pulsing(1,end,icell) = sum(egfmeki.dists(egfmeki.celltype == isite) > puls_thres) / sum(egfmeki.celltype == isite);
        egfmeki_pulsing_bd(:,end,icell) = bootci(2000,{ratio_fun,egfmeki.dists(egfmeki.celltype == isite)},'alpha',.32);
    end
    
end


for i = 1:length(sites_igfakti)
    
    isite = sites_igfakti(i);
    s = siteprop(isite,extension_igfakti);
    
    icell = 0;
    if ~isempty(strmatch(s.celltype,'MCF10A','exact'))
        icell = 1;
    elseif ~isempty(strmatch(s.celltype,'184A1','exact'))
        icell = 2;
    end
    
    if s.lig_dose == igf_dose && icell > 0
        igfakti_pulsing(1,s.drug_dose == akti_doses,icell) = sum(igfakti.dists(igfakti.celltype == isite) > puls_thres) / sum(igfakti.celltype == isite);
        igfakti_pulsing_bd(:,s.drug_dose == akti_doses,icell) = bootci(2000,{ratio_fun,igfakti.dists(igfakti.celltype == isite)},'alpha',.32);
    elseif s.lig_dose == 0 && s.drug_dose == 0 && icell > 0  
        igfakti_pulsing(1,end,icell) = sum(igfakti.dists(igfakti.celltype == isite) > puls_thres) / sum(igfakti.celltype == isite);
        igfakti_pulsing_bd(:,end,icell) = bootci(2000,{ratio_fun,igfakti.dists(igfakti.celltype == isite)},'alpha',.32);
    end
    
end

cell_names = {'MCF10A','184A1'};
ylim = [0 .4; 0 .8];
for icell = 1:2

    figure
    subplot(1,2,1)
    bar(1:length(meki_doses)+1,fliplr(egfmeki_pulsing(1,:,icell)),'b')
    hold on
    h = errorbar(1:length(meki_doses)+1,fliplr(egfmeki_pulsing(1,:,icell)),fliplr(egfmeki_pulsing(1,:,icell)-egfmeki_pulsing_bd(1,:,icell)),fliplr(egfmeki_pulsing_bd(2,:,icell)-egfmeki_pulsing(1,:,icell)),'k'); set(h,'linestyle','none')
    set(gca,'YLim',ylim(icell,:))
    set(gca,'XLim',[.5 length(meki_doses)+1.5])
    set(gca,'XTick',1:length(meki_doses)+1,'XTickLabel',[{'NS'},num2cell(fliplr(meki_doses))])
    xlabel('MEKi [muM]')
    ylabel('Fraction of pulsing cells')
    title(sprintf('%s: EGF %i ng/mL',cell_names{icell},egf_dose))
    subplot(1,2,2)
    bar(1:length(akti_doses)+1,fliplr(igfakti_pulsing(1,:,icell)),'r')
    hold on
    h = errorbar(1:length(akti_doses)+1,fliplr(igfakti_pulsing(1,:,icell)),fliplr(igfakti_pulsing(1,:,icell)-igfakti_pulsing_bd(1,:,icell)),fliplr(igfakti_pulsing_bd(2,:,icell)-igfakti_pulsing(1,:,icell)),'k'); set(h,'linestyle','none')
    set(gca,'YLim',ylim(icell,:))
    set(gca,'XLim',[.5 length(akti_doses)+1.5])
    xlabel('AKTi [muM]')
    set(gca,'XTick',1:length(akti_doses)+1,'XTickLabel',[{'NS'},num2cell(fliplr(akti_doses))])
    title(sprintf('%s: IGF %i ng/mL',cell_names{icell},igf_dose))
    
end