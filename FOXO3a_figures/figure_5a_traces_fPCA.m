% Figure 5a: Traces and fPCA for varying AKTi / MEKi dose
onlyellipses = 0; % Only plot ellipses for colored conditions; If onlyellipses = 0 also dots are plotted
markersizecolored = 6; % Specify size of colored markers
markersizegray = 3; % Specify size of gray background markers

egf_dose = 20; % Choose EGF dose here
igf_dose = 100; % Same for IGF

smooth_level = 500; % the larger the more smoothing

addpath('./Functions/')
egfmeki = load('Workspaces/dists_04182014');
egfmeki_fPCA = load('Workspaces/scores_04182014');
igfakti = load('./Workspaces/c_signal_03302014');

sites_egfmeki_unsorted = [1:39 41:72];
sites_egfmeki = [6:-1:1 7:12 19:24 18:-1:13 30:-1:25 31:39 41:72];
sites_igfakti_unsorted = 1:72;
sites_igfakti = sites_igfakti_unsorted([60:-1:1 61:72]);

extension_egfmeki = '04-18-2014';
extension_igfakti = '03-30-2014';

meki_doses = [.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0];
akti_doses = [.5 .1 .025 .00625 0];

pcs = [2 3];

figure

colmap = jet(length(meki_doses));
cell_names = {'MCF10A','184A1'};
legh = [];
legstr = {};
titstr = {};
for i = 1:length(sites_egfmeki)
    
    isite = sites_egfmeki(i);
    s = siteprop(isite,extension_egfmeki);
    
    icell = 0;
    if ~isempty(strmatch(s.celltype,'MCF10A','exact'))
        icell = 1;
    elseif ~isempty(strmatch(s.celltype,'184A1','exact'))
        icell = 2;
    end
    
    plotcond = 0;
    if ~isempty(strmatch(s.lig_name,'EGF','exact')) && s.lig_dose == egf_dose && icell > 0 
        mycol = colmap(meki_doses == s.drug_dose,:);
        plotcond = 1;
        mylegstr = sprintf('%s %g',s.drug_name, s.drug_dose);
        titstr{icell} = sprintf('%s: %s %i',cell_names{icell},s.lig_name,s.lig_dose);
    elseif s.lig_dose == 0 && s.drug_dose == 0 && icell > 0
        mycol = [.4 .4 .4];
        plotcond = 1;
        mylegstr = 'NS';
    end
    
    if plotcond
        subplot(1,2,icell)
        hold on
        title(titstr{icell})
        set(gca,'XLim',[50 400],'XTick',120:100:320,'XTickLabel',0:100:200)
        set(gca,'YLim',[-.015 .025])
        
        mean_signal = csaps(egfmeki.timestamp,nanmean(egfmeki.c_signal(:,egfmeki.celltype == isite),2),1/smooth_level,egfmeki.timestamp);
        std_signal = csaps(egfmeki.timestamp,nanstd(egfmeki.c_signal(:,egfmeki.celltype == isite),[],2),1/smooth_level,egfmeki.timestamp);
        n_signal = sum(egfmeki.celltype == isite);
        
        if icell == 1
            legh = [legh plot(egfmeki.timestamp,mean_signal,'Color',mycol)];
            legstr{end+1} = mylegstr;
        else
            plot(egfmeki.timestamp,mean_signal,'Color',mycol)
        end
        
        tmpx = [egfmeki.timestamp; flipud(egfmeki.timestamp)];
        tmpy = [mean_signal + 2*std_signal/sqrt(n_signal); flipud(mean_signal - 2*std_signal/sqrt(n_signal))];

        ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
        set(ltmp, 'FaceColor', mycol*0.1+0.9, 'EdgeColor', mycol*0.1+0.9);
        ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
        set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', mycol*0.3+0.7);
    end
    
end

subplot(1,2,1)
legend(legh,legstr)
xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');


figure

colmap = jet(length(akti_doses));
cell_names = {'MCF10A','184A1'};
legh = [];
legstr = {};
titstr = {};
for i = 1:length(sites_igfakti)
    
    isite = sites_igfakti(i);
    s = siteprop(isite,extension_igfakti);
    
    icell = 0;
    if ~isempty(strmatch(s.celltype,'MCF10A','exact'))
        icell = 1;
    elseif ~isempty(strmatch(s.celltype,'184A1','exact'))
        icell = 2;
    end
    
    plotcond = 0;
    if ~isempty(strmatch(s.lig_name,'IGF','exact')) && s.lig_dose == igf_dose && icell > 0 
        mycol = colmap(akti_doses == s.drug_dose,:);
        plotcond = 1;
        mylegstr = sprintf('%s %g',s.drug_name, s.drug_dose);
        titstr{icell} = sprintf('%s: %s %i',cell_names{icell},s.lig_name,s.lig_dose);
    elseif s.lig_dose == 0 && s.drug_dose == 0 && icell > 0
        mycol = [.4 .4 .4];
        plotcond = 1;
        mylegstr = 'NS';
    end
    
    if plotcond
        subplot(1,2,icell)
        hold on
        title(titstr{icell})
        set(gca,'XLim',[50 400],'XTick',120:100:320,'XTickLabel',0:100:200)
        set(gca,'YLim',[-.025 .035])
        
        mean_signal = csaps(igfakti.timestamp,nanmean(igfakti.c_signal(:,igfakti.celltype == isite),2),1/smooth_level,igfakti.timestamp);
        std_signal = csaps(igfakti.timestamp,nanstd(igfakti.c_signal(:,igfakti.celltype == isite),[],2),1/smooth_level,igfakti.timestamp);
        n_signal = sum(igfakti.celltype == isite);
        
        if icell == 1
            legh = [legh plot(igfakti.timestamp,mean_signal,'Color',mycol)];
            legstr{end+1} = mylegstr;
        else
            plot(igfakti.timestamp,mean_signal,'Color',mycol)
        end
        
        tmpx = [igfakti.timestamp; flipud(igfakti.timestamp)];
        tmpy = [mean_signal + 2*std_signal/sqrt(n_signal); flipud(mean_signal - 2*std_signal/sqrt(n_signal))];

        ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
        set(ltmp, 'FaceColor', mycol*0.1+0.9, 'EdgeColor', mycol*0.1+0.9);
        ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
        set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', mycol*0.3+0.7);
    end
    
end

subplot(1,2,1)
legend(legh,legstr)
xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');


figure

for isb = 1:2
    subplot(1,2,isb)
    plot(egfmeki_fPCA.scores_all(pcs(1),:),egfmeki_fPCA.scores_all(pcs(2),:),'o','Color',[.7 .7 .7],'MarkerFaceColor',[.7 .7 .7],'MarkerSize',markersizegray)
    set(gca,'XLim',[-.25 .25],'YLim',[-.1 .1])
end

colmap = jet(length(meki_doses));
cell_names = {'MCF10A','184A1'};
legh = [];
legstr = {};
titstr = {};
for i = 1:length(sites_egfmeki)
    
    isite = sites_egfmeki(i);
    s = siteprop(isite,extension_egfmeki);
    
    icell = 0;
    if ~isempty(strmatch(s.celltype,'MCF10A','exact'))
        icell = 1;
    elseif ~isempty(strmatch(s.celltype,'184A1','exact'))
        icell = 2;
    end
    
    plotcond = 0;
    if ~isempty(strmatch(s.lig_name,'EGF','exact')) && s.lig_dose == egf_dose && icell > 0 
        mycol = colmap(meki_doses == s.drug_dose,:);
        plotcond = 1;
        mylegstr = sprintf('%s %g',s.drug_name, s.drug_dose);
        titstr{icell} = sprintf('%s: %s %i',cell_names{icell},s.lig_name,s.lig_dose);
    elseif s.lig_dose == 0 && s.drug_dose == 0 && icell > 0
        mycol = [.4 .4 .4];
        plotcond = 1;
        mylegstr = 'NS';
    end
    
    if plotcond
        subplot(1,2,icell)
        hold on
        title(titstr{icell})
        if ~onlyellipses
            plot(egfmeki_fPCA.scores_all(pcs(1),egfmeki_fPCA.celltype == isite),egfmeki_fPCA.scores_all(pcs(2),egfmeki_fPCA.celltype == isite),'o','Color',mycol,'MarkerFaceColor',mycol,'MarkerSize',markersizecolored)
        end
        plotEllipsis(egfmeki_fPCA.scores_all(pcs(1),egfmeki_fPCA.celltype == isite),egfmeki_fPCA.scores_all(pcs(2),egfmeki_fPCA.celltype == isite),mycol,.5);
    end
end

subplot(1,2,1)
xlabel(sprintf('PC %i',pcs(1)))
ylabel(sprintf('PC %i',pcs(2)))


figure

for isb = 1:2
    subplot(1,2,isb)
    plot(igfakti.scores_all(pcs(1),:),igfakti.scores_all(pcs(2),:),'o','Color',[.7 .7 .7],'MarkerFaceColor',[.7 .7 .7],'MarkerSize',markersizegray)
    set(gca,'XLim',[-.8 .5],'YLim',[-.2 .2])
end

colmap = jet(length(akti_doses));
cell_names = {'MCF10A','184A1'};
legh = [];
legstr = {};
titstr = {};
for i = 1:length(sites_igfakti)
    
    isite = sites_igfakti(i);
    s = siteprop(isite,extension_igfakti);
    
    icell = 0;
    if ~isempty(strmatch(s.celltype,'MCF10A','exact'))
        icell = 1;
    elseif ~isempty(strmatch(s.celltype,'184A1','exact'))
        icell = 2;
    end
    
    plotcond = 0;
    if ~isempty(strmatch(s.lig_name,'IGF','exact')) && s.lig_dose == igf_dose && icell > 0 
        mycol = colmap(akti_doses == s.drug_dose,:);
        plotcond = 1;
        mylegstr = sprintf('%s %g',s.drug_name, s.drug_dose);
        titstr{icell} = sprintf('%s: %s %i',cell_names{icell},s.lig_name,s.lig_dose);
    elseif s.lig_dose == 0 && s.drug_dose == 0 && icell > 0
        mycol = [.4 .4 .4];
        plotcond = 1;
        mylegstr = 'NS';
    end
    
    if plotcond
        subplot(1,2,icell)
        hold on
        title(titstr{icell})
        if ~onlyellipses
            plot(igfakti.scores_all(pcs(1),igfakti.celltype == isite),igfakti.scores_all(pcs(2),igfakti.celltype == isite),'o','Color',mycol,'MarkerFaceColor',mycol,'MarkerSize',markersizecolored)
        end
        plotEllipsis(igfakti.scores_all(pcs(1),igfakti.celltype == isite),igfakti.scores_all(pcs(2),igfakti.celltype == isite),mycol,.5);
    end
end

subplot(1,2,1)
xlabel(sprintf('PC %i',pcs(1)))
ylabel(sprintf('PC %i',pcs(2)))