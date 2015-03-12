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

puls_thres = [.2 .5; .3 .3];

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

zeromeki = find(meki_doses == 0);
zeroakti = find(akti_doses == 0);

egfmeki_pvals = nan*egfmeki_pulsing;
igfakti_pvals = nan*igfakti_pulsing;

for iligdose = 1:size(egfmeki_pulsing,1)
    for icell = 1:2
        refcond_mean = egfmeki_pulsing(iligdose,zeromeki,icell); % MEAN2
        refcond_std = egfmeki_pulsing_bd(2*iligdose,zeromeki,icell) - egfmeki_pulsing_bd(2*iligdose-1,zeromeki,icell); % STD2
        for idrugdose = 1:size(egfmeki_pulsing,2)
            if idrugdose ~= zeromeki
%                 egfmeki_pulsing(iligdose,idrugdose,icell); % MEAN1
%                 egfmeki_pulsing_bd(2*iligdose,idrugdose,icell) - egfmeki_pulsing_bd(2*iligdose-1,idrugdose,icell); % STD1

                egfmeki_pvals(iligdose,idrugdose,icell) = 1-cdf('t',abs((refcond_mean - egfmeki_pulsing(iligdose,idrugdose,icell)) / sqrt((refcond_std)^2 + (egfmeki_pulsing_bd(2*iligdose,idrugdose,icell) - egfmeki_pulsing_bd(2*iligdose-1,idrugdose,icell))^2)),1e6);
            end
        end
    end
end

for iligdose = 1:size(igfakti_pulsing,1)
    for icell = 1:2
        refcond_mean = igfakti_pulsing(iligdose,zeroakti,icell); % MEAN2
        refcond_std = igfakti_pulsing_bd(2*iligdose,zeroakti,icell) - igfakti_pulsing_bd(2*iligdose-1,zeroakti,icell); % STD2
        for idrugdose = 1:size(igfakti_pulsing,2)
            if idrugdose ~= zeroakti
%                 igfakti_pulsing(iligdose,idrugdose,icell); % MEAN1
%                 igfakti_pulsing_bd(2*iligdose,idrugdose,icell) - igfakti_pulsing_bd(2*iligdose-1,idrugdose,icell); % STD1

                igfakti_pvals(iligdose,idrugdose,icell) = 1-cdf('t',abs((refcond_mean - igfakti_pulsing(iligdose,idrugdose,icell)) / sqrt((refcond_std)^2 + (igfakti_pulsing_bd(2*iligdose,idrugdose,icell) - igfakti_pulsing_bd(2*iligdose-1,idrugdose,icell))^2)),1e6);
            end
        end
    end
end

egfmeki_signif1 = egfmeki_pvals < .1;
egfmeki_signif1 = double(egfmeki_signif1);
egfmeki_signif1(egfmeki_signif1 == 0) = nan;
egfmeki_signif2 = egfmeki_pvals < .05;
egfmeki_signif2 = double(egfmeki_signif2);
egfmeki_signif2(egfmeki_signif2 == 0) = nan;
igfakti_signif1 = igfakti_pvals < .01;
igfakti_signif1 = double(igfakti_signif1);
igfakti_signif1(igfakti_signif1 == 0) = nan;
igfakti_signif2 = igfakti_pvals < .05;
igfakti_signif2 = double(igfakti_signif2);
igfakti_signif2(igfakti_signif2 == 0) = nan;
    

cell_names = {'MCF10A','184A1'};
ylim = [0 1; 0 1];
matVer = ver('MATLAB');
matVer = str2double(matVer.Version)>=8.4;
for icell = 1:2

    figure
    subplot(1,2,1)
    h = bar(fliplr(egfmeki_pulsing(:,:,icell))');
    set(h(1),'FaceColor',[.3 .3 1])
    set(h(2),'FaceColor',[.15 .15 .5])
    hold on
    for ih = 1:length(h)
        if ~matVer
            x1 = get(get(h(ih),'children'), 'xdata');
            x1 = mean(x1([1 3],:));
        else
            x1 = bsxfun(@plus, h(ih).XData, [h(ih).XOffset]');
        end
        h1 = errorbar(x1,fliplr(egfmeki_pulsing(ih,:,icell)),fliplr(egfmeki_pulsing(ih,:,icell)-egfmeki_pulsing_bd(1+2*(ih-1),:,icell)),fliplr(egfmeki_pulsing_bd(2*ih,:,icell)-egfmeki_pulsing(ih,:,icell)),'k');
        set(h1,'linestyle','none')
        plot(x1,fliplr(egfmeki_signif1(ih,:,icell)) * .85,'k*')
        plot(x1,fliplr(egfmeki_signif2(ih,:,icell)) * .87,'k*')
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
    set(h(1),'FaceColor',[1 .3 .3])
    set(h(2),'FaceColor',[.5 .15 .15])
    hold on
    for ih = 1:length(h)
        if ~matVer
            x1 = get(get(h(ih),'children'), 'xdata');
            x1 = mean(x1([1 3],:));
        else
            x1 = bsxfun(@plus, h(ih).XData, [h(ih).XOffset]');
        end
        h1 = errorbar(x1,fliplr(igfakti_pulsing(ih,:,icell)),fliplr(igfakti_pulsing(ih,:,icell)-igfakti_pulsing_bd(1+2*(ih-1),:,icell)),fliplr(igfakti_pulsing_bd(2*ih,:,icell)-igfakti_pulsing(ih,:,icell)),'k');
        set(h1,'linestyle','none')
        plot(x1,fliplr(igfakti_signif1(ih,:,icell)) * .85,'k*')
        plot(x1,fliplr(igfakti_signif2(ih,:,icell)) * .87,'k*')
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
