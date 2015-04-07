% Figure 5c: Bar charts EGF4 varying MEKi and IGF varying AKTi
egf_dose = [4 20]; % Choose EGF dose here
igf_dose = [20 100]; % Same for IGF

addpath('./Functions/')
egfmeki = load('Workspaces/scores_04182014');
igfakti = load('./Workspaces/c_signal_03302014');

pc = 2;

sites_egfmeki = [1:39 41:72];
sites_igfakti = 1:72;

extension_egfmeki = '04-18-2014';
extension_igfakti = '03-30-2014';

puls_thres = [.2 .5; .3 .3];

meki_doses = [.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0];
akti_doses = [.5 .1 .025 .00625 0];

egfmeki_pulsing = nan(length(egf_dose),length(meki_doses)+1,2);
igfakti_pulsing = nan(length(igf_dose),length(akti_doses)+1,2);

egfmeki_pulsing_bd = nan(length(egf_dose),length(meki_doses)+1,2);
igfakti_pulsing_bd = nan(length(igf_dose),length(akti_doses)+1,2);

rng(0) % Make sure that bootstrap samples are reproducible
ratio_fun = @(x) nanmean(x);

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
            egfmeki_pulsing(id,s.drug_dose == meki_doses,icell) = nanmean(egfmeki.scores_all(pc,egfmeki.celltype == isite));
            egfmeki_pulsing_bd(id,s.drug_dose == meki_doses,icell) = nanstd(egfmeki.scores_all(pc,egfmeki.celltype == isite))./sqrt(sum(egfmeki.celltype == isite));
        elseif s.lig_dose == 0 && s.drug_dose == 0 && icell > 0  
            egfmeki_pulsing(id,end,icell) = nanmean(egfmeki.scores_all(pc,egfmeki.celltype == isite));
            egfmeki_pulsing_bd(id,end,icell) = nanstd(egfmeki.scores_all(pc,egfmeki.celltype == isite))./sqrt(sum(egfmeki.celltype == isite));
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
            igfakti_pulsing(id,s.drug_dose == akti_doses,icell) = nanmean(igfakti.scores_all(pc,igfakti.celltype == isite));
            igfakti_pulsing_bd(id,s.drug_dose == akti_doses,icell) = nanstd(igfakti.scores_all(pc,igfakti.celltype == isite))./sqrt(sum(igfakti.celltype == isite));
        elseif s.lig_dose == 0 && s.drug_dose == 0 && icell > 0  
            igfakti_pulsing(id,end,icell) = nanmean(igfakti.scores_all(pc,igfakti.celltype == isite));
            igfakti_pulsing_bd(id,end,icell) = nanstd(igfakti.scores_all(pc,igfakti.celltype == isite))./sqrt(sum(igfakti.celltype == isite));
        end

    end
    
end

cell_names = {'MCF10A','184A1'};
ylim = [-.12 .25; -.12 .25];
matVer = ver('MATLAB');
matVer = str2double(matVer.Version)>=8.4;
for icell = 1:2

    figure
    subplot(1,2,1)
    h = bar(fliplr(egfmeki_pulsing(1,[1 6 7],icell))');
    set(h(1),'FaceColor',[.3 .3 1])
%     set(h(2),'FaceColor',[.15 .15 .5])
    hold on
    for ih = 1:length(h)
        if ~matVer
            x1 = get(get(h(ih),'children'), 'xdata');
            x1 = mean(x1([1 3],:));
        else
            x1 = bsxfun(@plus, h(ih).XData, [h(ih).XOffset]');
        end
        h1 = errorbar(x1,fliplr(egfmeki_pulsing(ih,[1 6 7],icell)),fliplr(egfmeki_pulsing_bd(ih,[1 6 7],icell)),'k');
        set(h1,'linestyle','none')
    end
    set(gca,'YLim',ylim(icell,:))
    set(gca,'XLim',[.5 2+1.5])
    set(gca,'XTick',1:3,'XTickLabel',[{'NS'},num2cell(fliplr(meki_doses([1 end])))])
    xlabel('MEKi [muM]')
    ylabel(sprintf('Mean PC%i score',pc))
    title(sprintf('%s: EGF 4 ng/mL',cell_names{icell}))
%     legstr = {};
%     for id = 1:length(egf_dose)
%         legstr{end+1} = sprintf('%i ng/mL',egf_dose(id));
%     end
%     legend(legstr)
    
    subplot(1,2,2)
    h = bar(fliplr(igfakti_pulsing(1,[1 5 6],icell))');
    set(h(1),'FaceColor',[1 .3 .3])
%     set(h(2),'FaceColor',[.5 .15 .15])
    hold on
    for ih = 1:length(h)
        if ~matVer
            x1 = get(get(h(ih),'children'), 'xdata');
            x1 = mean(x1([1 3],:));
        else
            x1 = bsxfun(@plus, h(ih).XData, [h(ih).XOffset]');
        end
        h1 = errorbar(x1,fliplr(igfakti_pulsing(ih,[1 5 6],icell)),fliplr(igfakti_pulsing_bd(ih,[1 5 6],icell)),'k');
        set(h1,'linestyle','none')
    end
    set(gca,'YLim',ylim(icell,:))
    set(gca,'XLim',[.5 2+1.5])
    xlabel('AKTi [muM]')
    set(gca,'XTick',1:3,'XTickLabel',[{'NS'},num2cell(fliplr(akti_doses([1 end])))])
    title(sprintf('%s: IGF 20 ng/mL',cell_names{icell}))
%     legstr = {};
%     for id = 1:length(igf_dose)
%         legstr{end+1} = sprintf('%i ng/mL',igf_dose(id));
%     end
%     legend(legstr)
    
end
