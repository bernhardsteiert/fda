% Figure 5d: Pulsing vs. early PC2
egf_dose = [4 20]; % Choose EGF dose here
igf_dose = [20 100]; % Same for IGF

pc = 2;

addpath('./Functions/')
egfmeki = load('Workspaces/dists_04182014');
egfmeki_fPCA = load('Workspaces/scores_04182014');
igfakti = load('./Workspaces/c_signal_03302014');

extension_egfmeki = '04-18-2014';
extension_igfakti = '03-30-2014';

sites_egfmeki_unsorted = [1:39 41:72];
sites_egfmeki = [6:-1:1 7:12 19:24 18:-1:13 30:-1:25 31:39 41:72];
sites_igfakti_unsorted = 1:72;
sites_igfakti = sites_igfakti_unsorted([60:-1:1 61:72]);

puls_thres = [.2 .5; .3 .3];

meki_doses = [.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0];
akti_doses = [.5 .1 .025 .00625 0];

egfmeki_pulsing = nan(length(egf_dose),length(meki_doses)+1,2);
igfakti_pulsing = nan(length(igf_dose),length(akti_doses)+1,2);
egfmeki_early = egfmeki_pulsing;
igfakti_early = igfakti_pulsing;

egfmeki_pulsing_gray = cell(1,2);
egfmeki_early_gray = cell(1,2);
igfakti_pulsing_gray = cell(1,2);
igfakti_early_gray = cell(1,2);


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
            egfmeki_early(id,s.drug_dose == meki_doses,icell) = nanmean(egfmeki_fPCA.scores_all(pc,egfmeki.celltype == isite));
%         elseif s.lig_dose == 0 && s.drug_dose == 0 && icell > 0  
%             egfmeki_pulsing(id,end,icell) = sum(egfmeki.dists(egfmeki.celltype == isite) > puls_thres(icell,1)) / sum(egfmeki.celltype == isite);
%             egfmeki_early(id,end,icell) = nanmean(egfmeki_fPCA.scores_all(pc,egfmeki.celltype == isite));
        elseif icell > 0 && ~isempty(strmatch(s.lig_name,'EGF','exact'))
            egfmeki_pulsing_gray{icell} = [egfmeki_pulsing_gray{icell} sum(egfmeki.dists(egfmeki.celltype == isite) > puls_thres(icell,1)) / sum(egfmeki.celltype == isite)];
            egfmeki_early_gray{icell} = [egfmeki_early_gray{icell} nanmean(egfmeki_fPCA.scores_all(pc,egfmeki.celltype == isite))];
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
            igfakti_early(id,s.drug_dose == akti_doses,icell) = nanmean(igfakti.scores_all(pc,igfakti.celltype == isite));
%         elseif s.lig_dose == 0 && s.drug_dose == 0 && icell > 0  
%             igfakti_pulsing(id,end,icell) = sum(igfakti.dists(igfakti.celltype == isite) > puls_thres(icell,1)) / sum(igfakti.celltype == isite);
%             igfakti_early(id,end,icell) = nanmean(igfakti_fPCA.scores_all(pc,igfakti.celltype == isite));
        elseif icell > 0 && ~isempty(strmatch(s.lig_name,'IGF','exact'))
            igfakti_pulsing_gray{icell} = [igfakti_pulsing_gray{icell} sum(igfakti.dists(igfakti.celltype == isite) > puls_thres(icell,2)) / sum(igfakti.celltype == isite)];
            igfakti_early_gray{icell} = [igfakti_early_gray{icell} nanmean(igfakti.scores_all(pc,igfakti.celltype == isite))];
        end

    end
    
end



cell_names = {'MCF10A','184A1'};
for icell = 1:2
    figure
    subplot(1,3,1)
    
    plot(egfmeki_early_gray{icell},egfmeki_pulsing_gray{icell},'^','Color',[.7 .7 .7]);
    hold on
    colmap = jet(length(meki_doses));
    colmap = [colmap; 0 0 0];
    markers = {'s','o'};
    for i = 1:size(egfmeki_early,2)
        col = colmap(i,:);
        legh = [];
        legstr = {};
        for j = 1:size(egfmeki_early,1)
            legh = [legh plot(egfmeki_early(j,i,icell),egfmeki_pulsing(j,i,icell),markers{j},'MarkerFaceColor',col)];
            legstr{end+1} = sprintf('EGF %g ng/mL',egf_dose(j));
        end
    end
    legend(legh,legstr)
    title(sprintf('%s: %s / MEKi',cell_names{icell},'EGF'))
    set(gca,'XLim',[-.2 .3],'YLim',[0 1])
    xlabel(sprintf('Early PC%i',pc))
    ylabel('Fraction of pulsing cells')
    
    subplot(1,3,2)
    
    plot(igfakti_early_gray{icell},igfakti_pulsing_gray{icell},'^','Color',[.7 .7 .7]);
    hold on
    colmap = jet(length(akti_doses));
    colmap = [colmap; 0 0 0];
    markers = {'s','o'};
    for i = 1:size(igfakti_early,2)
        col = colmap(i,:);
        legh = [];
        legstr = {};
        for j = 1:size(igfakti_early,1)
            legh = [legh plot(igfakti_early(j,i,icell),igfakti_pulsing(j,i,icell),markers{j},'MarkerFaceColor',col)];
            legstr{end+1} = sprintf('IGF %g ng/mL',igf_dose(j));
        end
    end
    legend(legh,legstr)
    title(sprintf('%s: %s / AKTi',cell_names{icell},'IGF'))
    set(gca,'XLim',[-.2 .3],'YLim',[0 1])
    xlabel(sprintf('Early PC%i',pc))
    ylabel('Fraction of pulsing cells')
    
    subplot(1,3,3)
    set(gca,'Visible','off','CLim',[0 1])
    colormap(jet)
    c = colorbar('Location','West');
    set(c,'YTick',[0 .5 1],'YTickLabel',{'High drug dose','Medium drug dose','No Drug'});
    
end