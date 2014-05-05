% Figure 5b: EGF vs. MEKi Traces
addpath('./Functions/')
load('Workspaces/dists_04182014')
extension = '04-18-2014';
sites_all = [1:39 41:72];

figure

smooth_level = 500; % the larger the more smoothing
mean_signal = nan(length(timestamp),sites_all(end));
std_signal = nan(length(timestamp),sites_all(end));
n_signal = nan(length(timestamp),sites_all(end));
for isite = sites_all
    mean_signal(:,isite) = csaps(timestamp,nanmean(c_signal(:,celltype == isite),2),1/smooth_level,timestamp);
    std_signal(:,isite) = csaps(timestamp,nanstd(c_signal(:,celltype == isite),[],2),1/smooth_level,timestamp);
    n_signal(isite) = length(c_signal(:,celltype == isite));
end

fixedEGF = 13:24;
fixedMEKi = [4 24-4+1 24+4 48-4+1 48+4 72-4+1 9 24-9+1 24+9 48-9+1 48+9 72-9+1];
fixedEGF = fixedEGF([6:-1:1 7:12]);
fixedMEKi = fixedMEKi([6:-1:1 7:12]);

ligand_dose = [100 20 4 .8 .16 0];
drug_dose = [.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0];

colmap = jet(length(drug_dose));
linewidth = 1;
legh = [];
legstr = {};
for isite = fixedEGF
    sprop = siteprop(isite,extension);
    if strmatch(sprop.celltype,'MCF10A','exact')
        subplot(2,2,1)
        hold on
        legh = [legh plot(timestamp,mean_signal(:,isite),'Color',colmap(sprop.drug_dose == drug_dose,:))];
        legstr{end+1} = sprintf('MEKi %g',drug_dose(sprop.drug_dose == drug_dose));
    else
        subplot(2,2,2)
        hold on
        plot(timestamp,mean_signal(:,isite),'Color',colmap(sprop.drug_dose == drug_dose,:))
    end
    
    tmpx = [timestamp; flipud(timestamp)];
    tmpy = [mean_signal(:,isite) + 2*std_signal(:,isite)/sqrt(n_signal(isite)); flipud(mean_signal(:,isite) - 2*std_signal(:,isite)/sqrt(n_signal(isite)))];
    
    ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp, 'FaceColor', colmap(sprop.drug_dose == drug_dose,:)*0.1+0.9, 'EdgeColor', colmap(sprop.drug_dose == drug_dose,:)*0.1+0.9,'LineWidth',linewidth);
    ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', colmap(sprop.drug_dose == drug_dose,:)*0.3+0.7,'LineWidth',linewidth);
    
end
legend(legh,legstr)

legh = [];
legstr = {};
colmap = jet(length(ligand_dose));
for isite = fixedMEKi
    sprop = siteprop(isite,extension);
    if strmatch(sprop.celltype,'MCF10A','exact')
        subplot(2,2,3)
        hold on
        legh = [legh plot(timestamp,mean_signal(:,isite),'Color',colmap(sprop.lig_dose == ligand_dose,:))];
        legstr{end+1} = sprintf('EGF %g',ligand_dose(sprop.lig_dose == ligand_dose));
    else
        subplot(2,2,4)
        hold on
        plot(timestamp,mean_signal(:,isite),'Color',colmap(sprop.lig_dose == ligand_dose,:))
    end
    
    tmpx = [timestamp; flipud(timestamp)];
    tmpy = [mean_signal(:,isite) + 2*std_signal(:,isite)/sqrt(n_signal(isite)); flipud(mean_signal(:,isite) - 2*std_signal(:,isite)/sqrt(n_signal(isite)))];
    
    ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp, 'FaceColor', colmap(sprop.lig_dose == ligand_dose,:)*0.1+0.9, 'EdgeColor', colmap(sprop.lig_dose == ligand_dose,:)*0.1+0.9,'LineWidth',linewidth);
    ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', colmap(sprop.lig_dose == ligand_dose,:)*0.3+0.7,'LineWidth',linewidth);
    
end
legend(legh,legstr)

subplot(2,2,1)
title('MCF10A - EGF 20; varying MEKi')
set(gca,'YLim',[-.012 .022])
set(gca,'XLim',[0 1000])
ylabel('log_{10} FOXO3a [Cyt/Nuc]');
subplot(2,2,2)
title('184A1 - EGF 20; varying MEKi')
set(gca,'YLim',[-.012 .022])
set(gca,'XLim',[0 1000])
subplot(2,2,3)
title('MCF10A - MEKi 1/640; varying EGF')
set(gca,'YLim',[-.012 .022])
set(gca,'XLim',[0 1000])
xlabel('time [min]')
subplot(2,2,4)
title('184A1 - MEKi 1/640; varying EGF')
set(gca,'YLim',[-.012 .022])
set(gca,'XLim',[0 1000])