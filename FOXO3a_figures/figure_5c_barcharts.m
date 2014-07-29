% Figure 5c: Bar charts EGF4 varying MEKi and IGF varying AKTi
egf_dose = [4 20]; % Choose EGF dose here
igf_dose = [20 100]; % Same for IGF

addpath('./Functions/')
egfmeki = load('Workspaces/dists_04182014');
igfakti = load('./Workspaces/c_signal_03302014');

sites_egfmeki = [1:39 41:72];
sites_igfakti = 1:72;

extension_egfmeki = '04-18-2014';
extension_igfakti = '03-30-2014';

puls_thres = [.3 .5; .3 .3];

meki_doses = [.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0];
akti_doses = [.5 .1 .025 .00625 0];

egfmeki_pulsing = nan(length(egf_dose),length(meki_doses)+1,2);
igfakti_pulsing = nan(length(igf_dose),length(akti_doses)+1,2);

egfmeki_pulsing_bd = nan(2*length(egf_dose),length(meki_doses)+1,2);
igfakti_pulsing_bd = nan(2*length(igf_dose),length(akti_doses)+1,2);

rng(0) % Make sure that bootstrap samples are reproducible
ratio_fun = @(x,puls_thres) sum(x > puls_thres) / length(x);

for id = 1:length(egf_dose)

    for i = 1:length(sites_egfmeki)

        isite = sites_egfmeki(i);
        s = siteprop(isite,extension_egfmeki);

        icell = 0;
        if ~isempty(strmatch(s.celltype,'MCF10A','exact'))
            icell = 1;
        elseif ~isempty(strmatch(s.celltype,'184A1','exact'))
            icell = 2;
        end

        if s.lig_dose == egf_dose(id) && icell > 0 && ~isempty(strmatch(s.lig_name,'EGF','exact'))
            egfmeki_pulsing(id,s.drug_dose == meki_doses,icell) = sum(egfmeki.dists(egfmeki.celltype == isite) > puls_thres(icell,1)) / sum(egfmeki.celltype == isite);
            egfmeki_pulsing_bd(1+2*(id-1):2*id,s.drug_dose == meki_doses,icell) = bootci(2000,{ratio_fun,egfmeki.dists(egfmeki.celltype == isite),puls_thres(icell,1)},'alpha',.32);
        elseif s.lig_dose == 0 && s.drug_dose == 0 && icell > 0  
            egfmeki_pulsing(id,end,icell) = sum(egfmeki.dists(egfmeki.celltype == isite) > puls_thres(icell,1)) / sum(egfmeki.celltype == isite);
            egfmeki_pulsing_bd(1+2*(id-1):2*id,end,icell) = bootci(2000,{ratio_fun,egfmeki.dists(egfmeki.celltype == isite),puls_thres(icell,1)},'alpha',.32);
        end

    end
    
end

for id = 1:length(igf_dose)

    for i = 1:length(sites_igfakti)

        isite = sites_igfakti(i);
        s = siteprop(isite,extension_igfakti);

        icell = 0;
        if ~isempty(strmatch(s.celltype,'MCF10A','exact'))
            icell = 1;
        elseif ~isempty(strmatch(s.celltype,'184A1','exact'))
            icell = 2;
        end

        if s.lig_dose == igf_dose(id) && icell > 0 && ~isempty(strmatch(s.lig_name,'IGF','exact'))
            igfakti_pulsing(id,s.drug_dose == akti_doses,icell) = sum(igfakti.dists(igfakti.celltype == isite) > puls_thres(icell,2)) / sum(igfakti.celltype == isite);
            igfakti_pulsing_bd(1+2*(id-1):2*id,s.drug_dose == akti_doses,icell) = bootci(2000,{ratio_fun,igfakti.dists(igfakti.celltype == isite),puls_thres(icell,2)},'alpha',.32);
        elseif s.lig_dose == 0 && s.drug_dose == 0 && icell > 0  
            igfakti_pulsing(id,end,icell) = sum(igfakti.dists(igfakti.celltype == isite) > puls_thres(icell,2)) / sum(igfakti.celltype == isite);
            igfakti_pulsing_bd(1+2*(id-1):2*id,end,icell) = bootci(2000,{ratio_fun,igfakti.dists(igfakti.celltype == isite),puls_thres(icell,2)},'alpha',.32);
        end

    end
    
end

cell_names = {'MCF10A','184A1'};
ylim = [0 1; 0 1];
for icell = 1:2

    figure
    subplot(1,2,1)
    h = bar(fliplr(egfmeki_pulsing(:,:,icell))');
    set(h(1),'FaceColor',[0 0 1])
    set(h(2),'FaceColor',[1 0 0])
    hold on
    for ih = 1:length(h)
        x1 = get(get(h(ih),'children'), 'xdata');
        x1 = mean(x1([1 3],:));
        h1 = errorbar(x1,fliplr(egfmeki_pulsing(ih,:,icell)),fliplr(egfmeki_pulsing(ih,:,icell)-egfmeki_pulsing_bd(1+2*(ih-1),:,icell)),fliplr(egfmeki_pulsing_bd(2*ih,:,icell)-egfmeki_pulsing(ih,:,icell)),'k');
        set(h1,'linestyle','none')
    end
    set(gca,'YLim',ylim(icell,:))
    set(gca,'XLim',[.5 length(meki_doses)+1.5])
    set(gca,'XTick',1:length(meki_doses)+1,'XTickLabel',[{'NS'},num2cell(fliplr(meki_doses))])
    xlabel('MEKi [muM]')
    ylabel('Fraction of pulsing cells')
    title(sprintf('%s: EGF',cell_names{icell}))
    legstr = {};
    for id = 1:length(egf_dose)
        legstr{end+1} = sprintf('%i ng/mL',egf_dose(id));
    end
    legend(legstr)
    
    subplot(1,2,2)
    h = bar(fliplr(igfakti_pulsing(:,:,icell))');
    set(h(1),'FaceColor',[0 0 1])
    set(h(2),'FaceColor',[1 0 0])
    hold on
    for ih = 1:length(h)
        x1 = get(get(h(ih),'children'), 'xdata');
        x1 = mean(x1([1 3],:));
        h1 = errorbar(x1,fliplr(igfakti_pulsing(ih,:,icell)),fliplr(igfakti_pulsing(ih,:,icell)-igfakti_pulsing_bd(1+2*(ih-1),:,icell)),fliplr(igfakti_pulsing_bd(2*ih,:,icell)-igfakti_pulsing(ih,:,icell)),'k');
        set(h1,'linestyle','none')
    end
    set(gca,'YLim',ylim(icell,:))
    set(gca,'XLim',[.5 length(akti_doses)+1.5])
    xlabel('AKTi [muM]')
    set(gca,'XTick',1:length(akti_doses)+1,'XTickLabel',[{'NS'},num2cell(fliplr(akti_doses))])
    title(sprintf('%s: IGF',cell_names{icell}))
    legstr = {};
    for id = 1:length(igf_dose)
        legstr{end+1} = sprintf('%i ng/mL',igf_dose(id));
    end
    legend(legstr)
    
end
