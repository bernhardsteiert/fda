% Set puls_thres = 0 if gating should be off
puls_thres = .5;
% puls_thres = 0;

sites_all = [4:10 17:-1:11 24:30 37:-1:31 44:50 57:-1:51 64:70];
names_all = {'EGF','IGF','FGF','HRG','HGF','EPR','BTC'};
reorder = [4 7 1 6 5 2 3];

% Comment the following in if EGF / BTC shall be replaced by MEKi variants
% sites_all = [1 5:10 17:-1:11 24:30 37:-1:31 44:50 57:-1:51 61 65:70]; % EGF / BTC with MEKi

doses = [0 2.5 5 10 20 50 100];
colmap = flipud(winter(length(doses)));

nrows = 7;
% ncols = 9; % nEdge; SNR; ...
ncols = 8;

medians = nan(length(doses),nrows,ncols);
singles = [];
features = {'Final score','nEdge','SNR','Amplitude','Pairwise','Peak duration','Peak distance','Ratio pulsing cells'};
single_features = {'Final score','nEdge','SNR','Amplitude','Pairwise','Peak duration','Peak distance','Ratio pulsing cells'};
celltype = [];

for isite = sites_all
    
    [radial_dists c_signal_tmp tmp2 nEdge SNR amp pw peakdur_mean peakdur_std peakdis_mean peakdis_std] = edge_snr_score_pw_distdur(isite);
    
    s = siteprop(isite);
    
    doseind = s.lig_dose == doses;
%     medians(doseind,reorder(s.lig_index),1) = median(radial_dists(radial_dists > puls_thres));
%     medians(doseind,reorder(s.lig_index),2) = median(nEdge(radial_dists > puls_thres));
%     medians(doseind,reorder(s.lig_index),3) = median(SNR(radial_dists > puls_thres));
%     medians(doseind,reorder(s.lig_index),4) = median(amp(radial_dists > puls_thres));
%     medians(doseind,reorder(s.lig_index),5) = median(pw(radial_dists > puls_thres));
%     medians(doseind,reorder(s.lig_index),6) = median(peakdur_mean(~isnan(peakdur_mean) & (radial_dists > puls_thres)));
%     medians(doseind,reorder(s.lig_index),7) = median(peakdis_mean(~isnan(peakdis_mean) & (radial_dists > puls_thres)));


    medians(doseind,reorder(s.lig_index),1) = median(radial_dists);
    medians(doseind,reorder(s.lig_index),2) = median(nEdge);
    medians(doseind,reorder(s.lig_index),3) = median(SNR);
    medians(doseind,reorder(s.lig_index),4) = median(amp);
    medians(doseind,reorder(s.lig_index),5) = median(pw);
    medians(doseind,reorder(s.lig_index),6) = median(peakdur_mean(~isnan(peakdur_mean)));
    medians(doseind,reorder(s.lig_index),7) = median(peakdis_mean(~isnan(peakdis_mean)));



    medians(doseind,reorder(s.lig_index),8) = sum(radial_dists > puls_thres)/length(radial_dists);
    
    peakdur_mean(isnan(peakdur_mean)) = nanmean(peakdur_mean);
    peakdis_mean(isnan(peakdis_mean)) = nanmean(peakdis_mean);
    
    % Gating for PCA (if differences of pulsing cells are investigated)
    radial_dists = radial_dists(radial_dists > puls_thres);
    nEdge = nEdge(radial_dists > puls_thres);
    SNR = SNR(radial_dists > puls_thres);
    amp = amp(radial_dists > puls_thres);
    pw = pw(radial_dists > puls_thres);
    peakdur_mean = peakdur_mean(radial_dists > puls_thres);
    peakdis_mean = peakdis_mean(radial_dists > puls_thres);
    celltype = [celltype; ones(sum(radial_dists > puls_thres),1)*isite];
    
    singles = [singles; [radial_dists' nEdge' SNR' amp' pw' peakdur_mean' peakdis_mean']];
    
end

%% Plotting (Regression)
close all

f1 = figure;

xfac = 1;
yfac = 1;
fontsize = 6;

setFigure(f1,xfac,yfac,fontsize)

for irow = 1:nrows
    for icol = 1:ncols
        
        subplot(nrows,ncols,(irow-1)*ncols+icol)
        hold on
        plot(medians(:,irow,icol),'x')
        
        [axb s] = polyfit(2:length(doses),medians(2:end,irow,icol)',1);
        plot(1:length(doses),(1:length(doses))*axb(1) + axb(2),'k--')
        
        if irow == 1
            title(features{icol})
        end
        if icol == 1
            ylabel(names_all{irow == reorder})
        end
        tmpmed = medians(:,:,icol);
        tmpmed = tmpmed(:);
        ylim = [min(tmpmed) max(tmpmed)]+[-1 1]*(range(tmpmed)+1e-10)*0.1;
        if ~isnan(ylim)
            set(gca,'YLim',ylim)
        end
        set(gca,'XLim',[.5 length(doses)+.5])
        set(gca,'XTick',1:length(doses),'XTickLabel',doses)
        
        plot([1.5 1.5],ylim,'r:')
        
        if irow == nrows && icol == 1
            xlabel('Ligand dose')
        end
    end
end

%% Plotting Regression only for ratio; Puls Dist; Amplitude
close all

f2 = figure;

xfac = 1;
yfac = .5;
fontsize = 10;

setFigure(f2,xfac,yfac,fontsize)

subplot(1,3,1)
hold on
icol = 8;
mycolor = lines(nrows);
legh = [];
for irow = 1:nrows
    plot(0:length(doses)-1,medians(:,irow,icol),'x','Color',mycolor(irow,:))

    [axb s] = polyfitZero(1:length(doses)-1,medians(2:end,irow,icol)'-mean(medians(1,:,icol)),1);
    legh = [legh plot(0:length(doses)-1,(0:length(doses)-1)*axb(1) + mean(medians(1,:,icol)),'-','Color',mycolor(irow,:))];
end
title(features{icol})
set(gca,'XLim',[-.5 length(doses)-.5])
set(gca,'XTick',0:length(doses)-1,'XTickLabel',doses)
set(gca,'YLim',[0 .5])
xlabel('Ligand dose')
legend(legh(reorder),names_all,'Location','NorthWest')

subplot(1,3,2)
hold on
icol = 4;
mycolor = lines(nrows);
legh = [];
for irow = 1:nrows
    plot(0:length(doses)-1,medians(:,irow,icol),'x','Color',mycolor(irow,:))

    [axb s] = polyfitZero(1:length(doses)-1,medians(2:end,irow,icol)'-mean(medians(1,:,icol)),1);
    legh = [legh plot(0:length(doses)-1,(0:length(doses)-1)*axb(1) + mean(medians(1,:,icol)),'-','Color',mycolor(irow,:))];
end
title(features{icol})
set(gca,'XLim',[-.5 length(doses)-.5])
set(gca,'XTick',0:length(doses)-1,'XTickLabel',doses)
% set(gca,'YLim',[7 10]*1e-3)
xlabel('Ligand dose')
% legend(legh(reorder),names_all,'Location','NorthWest')

subplot(1,3,3)
hold on
icol = 6;
mycolor = lines(nrows);
legh = [];
for irow = 1:nrows
    plot(0:length(doses)-1,medians(:,irow,icol),'x','Color',mycolor(irow,:))

    [axb s] = polyfit(1:length(doses)-1,medians(2:end,irow,icol)',1);
    legh = [legh plot(0:length(doses)-1,(0:length(doses)-1)*axb(1) + axb(2),'-','Color',mycolor(irow,:))];
end
title(features{icol})
set(gca,'XLim',[-.5 length(doses)-.5])
set(gca,'XTick',0:length(doses)-1,'XTickLabel',doses)
% set(gca,'YLim',[7 10]*1e-3)
xlabel('Ligand dose')
% legend(legh(reorder),names_all,'Location','NorthWest')

return


for irow = 1:nrows
    for icol = 1:ncols
        
        subplot(nrows,ncols,(irow-1)*ncols+icol)
        hold on
        plot(medians(:,irow,icol),'x')
        
        [axb s] = polyfit(2:length(doses),medians(2:end,irow,icol)',1);
        plot(1:length(doses),(1:length(doses))*axb(1) + axb(2),'k--')
        
        if irow == 1
            title(features{icol})
        end
        if icol == 1
            ylabel(names_all{irow == reorder})
        end
        tmpmed = medians(:,:,icol);
        tmpmed = tmpmed(:);
        ylim = [min(tmpmed) max(tmpmed)]+[-1 1]*(range(tmpmed)+1e-10)*0.1;
        if ~isnan(ylim)
            set(gca,'YLim',ylim)
        end
        set(gca,'XLim',[.5 length(doses)+.5])
        set(gca,'XTick',1:length(doses),'XTickLabel',doses)
        
        plot([1.5 1.5],ylim,'r:')
        
        if irow == nrows && icol == 1
            xlabel('Ligand dose')
        end
    end
end

%% Plotting (Densities)
close all

f3 = figure;

xfac = 1;
yfac = 1;
fontsize = 6;

setFigure(f3,xfac,yfac,fontsize)

for isite = sites_all
    s = siteprop(isite);
    irow = reorder(s.lig_index);
    for icol = 1:size(singles,2)
        
        subplot(nrows,ncols,(irow-1)*ncols+icol)
        hold on
        [f,xi] = ksdensity(singles(celltype == isite,icol),'width',range());
        plot(xi,f,'Color',colmap(s.lig_dose == doses,:))
        
        if irow == 1
            title(features{icol})
        end
        if icol == 1
            ylabel(names_all{irow == reorder})
        end
        tmpmed = medians(:,:,icol);
        tmpmed = tmpmed(:);
        
        if irow == nrows && icol == 1
            xlabel('Ligand dose')
        end
    end
end

%% Plotting (Correlations)
close all

f2 = figure;

xfac = 1;
yfac = 1;
fontsize = 6;

setFigure(f2,xfac,yfac,fontsize)

for irow = 1:size(singles,2)-1
    for icol = irow:size(singles,2)-1
        subplot(size(singles,2)-1,size(singles,2)-1,(irow-1)*(size(singles,2)-1)+icol)
        plot(singles(:,icol+1),singles(:,irow),'k.')
        hold on
        
        ylim = [min(singles(:,irow)) max(singles(:,irow))]+[-1 1]*range(singles(:,irow))*0.1;
        xlim = [min(singles(:,icol+1)) max(singles(:,icol+1))]+[-1 1]*range(singles(:,icol+1))*0.1;
        set(gca,'XLim',xlim,'YLim',ylim)
    end
end

for i = 1:size(singles,2)-1
    subplot(size(singles,2)-1,size(singles,2)-1,i)
    title(single_features{i+1})
    subplot(size(singles,2)-1,size(singles,2)-1,(size(singles,2))*(i-1)+1)
    ylabel(single_features{i})
end

%% Specifically investigate peak duration
close all

f2 = figure;

xfac = 1;
yfac = 1;
fontsize = 6;

setFigure(f2,xfac,yfac,fontsize)

for isite = sites_all
    s = siteprop(isite);
    irow = reorder(s.lig_index);
    icol = find(s.lig_dose == doses);
    subplot(length(names_all),length(doses),(irow-1)*length(doses) + icol)
    
    plot(singles(:,6),singles(:,1),'.','Color',[.7 .7 .7])
    hold on
    plot(singles(celltype == isite,6),singles(celltype == isite,1),'r.')
    
    ylabel('Pulsatory strength')
    xlabel(single_features{6})
    title(sprintf('%s %g',s.lig_name,s.lig_dose))
    
end


%% PCA for pulsatory behavior
close all

figure
hold on

inds = [2 3 4 6 7];
% inds = [2 3 4]; % wo Peak dur / Peak dist
[coeff scores latent] = princomp(singles(:,inds)./repmat(max(singles(:,inds),[],1),size(singles,1),1));

pcs = [2 3];

highdoses = [4 24 37 44 57 64 30];
% highdoses = [37 64 30];

color_ind = 1;
colmap = hsv(length(highdoses));
legstr = {};
for isite = highdoses
    s = siteprop(isite);
    titstr = s.lig_name;
    titstr = sprintf('%s %g',titstr,s.lig_dose);
    if s.inh_dose > 0
%         titstr = sprintf('%s%s %i muM',titstr,s.inh_name,s.inh_dose);
        titstr = sprintf('%s%s',titstr,s.inh_name);
    end
    legstr{end+1} = titstr;

    plot(scores(celltype == isite,pcs(1)),scores(celltype == isite,pcs(2)),'o','Color',colmap(isite == highdoses,:),'MarkerFaceColor',colmap(isite == highdoses,:))
    plotEllipsis(scores(celltype == isite,pcs(1)),scores(celltype == isite,pcs(2)),colmap(isite == highdoses,:),.5);
end

% xlim = [-.16 .24];
% set(gca,'XLim',xlim)

% axisEqual(get(gcf,'Position'))

ylabel(['PC ' num2str(pcs(2))])
% set(gca,'YTick',-.1:.1:.2)
xlabel(['PC ' num2str(pcs(1))])
% arrow([-.13 .06],[.02 .15],'Width',.5,'Length',7)

set(gca,'CLim',[0 1])
colormap(colmap)
colorbar('YTick',linspace(1./(2*length(highdoses)),1-1./(2*length(highdoses)),length(highdoses)),'YTickLabel',legstr,'TickLength', [0 0]) % Vertical colorbar
% return

%%
close all
figure

for ipc = 1:4
    
    subplot(2,2,ipc)
    plot(singles(:,1),scores(:,ipc),'k.')
    
    xlabel('Pulsatory strength')
    ylabel(['PC ' num2str(ipc)])
    
end

return

%% Linear Regression (depricated!)
close all

f3 = figure;

xfac = 1;
yfac = 1;
fontsize = 4;

setFigure(f3,xfac,yfac,fontsize)

for isite = sites_all
    
    [radial_dists c_signal_tmp tmp2 nEdge SNR amp pw peakdur_mean peakdur_std peakdis_mean peakdis_std] = edge_snr_score_pw_distdur(isite);
    
    s = siteprop(isite);
    
    if ~isempty(radial_dists)
        subplot(nrows,ncols,ncols*(reorder(s.lig_index)-1)+1)
        hold on
        plot(find(s.lig_dose == doses),median(radial_dists),'x')
        set(gca,'XLim',[.5 length(doses)+.5],'YLim',[-.01 .1])
        set(gca,'XTick',1:length(doses),'XTickLabel',doses)
        ylabel(s.lig_name)
        title('Final score')
    end
    
    if ~isempty(nEdge)
        subplot(nrows,ncols,ncols*(reorder(s.lig_index)-1)+2)
        hold on
        plot(find(s.lig_dose == doses),median(nEdge),'x')
        set(gca,'XLim',[.5 length(doses)+.5],'YLim',[-1 6])
        set(gca,'XTick',1:length(doses),'XTickLabel',doses)
        title('nEdge')
    end
    
    if ~isempty(SNR)
        subplot(nrows,ncols,ncols*(reorder(s.lig_index)-1)+3)
        hold on
        plot(find(s.lig_dose == doses),median(SNR),'x')
        set(gca,'XLim',[.5 length(doses)+.5],'YLim',[5 15])
        set(gca,'XTick',1:length(doses),'XTickLabel',doses)
        title('SNR')
    end
    
    if ~isempty(amp)
        subplot(nrows,ncols,ncols*(reorder(s.lig_index)-1)+4)
        hold on
        plot(find(s.lig_dose == doses),median(amp),'x')
        set(gca,'XLim',[.5 length(doses)+.5],'YLim',[0.002 0.01])
        set(gca,'XTick',1:length(doses),'XTickLabel',doses)
        title('Amplitude')
    end
    
    if ~isempty(pw)
        subplot(nrows,ncols,ncols*(reorder(s.lig_index)-1)+5)
        hold on
        plot(find(s.lig_dose == doses),median(pw),'x')
        set(gca,'XLim',[.5 length(doses)+.5],'YLim',[0 0.03])
        set(gca,'XTick',1:length(doses),'XTickLabel',doses)
        title('Pairwise')
    end
    
    if ~isempty(peakdur_mean)
        subplot(nrows,ncols,ncols*(reorder(s.lig_index)-1)+6)
        hold on
        plot(find(s.lig_dose == doses),median(peakdur_mean),'x')
        set(gca,'XLim',[.5 length(doses)+.5],'YLim',[50 160])
        set(gca,'XTick',1:length(doses),'XTickLabel',doses)
        title('Peak duration')
    end
    
%     if sum(~isnan(peakdur_std))
%         subplot(nrows,ncols,ncols*(reorder(s.lig_index)-1)+7)
%         hold on
%         [f,xi] = ksdensity(peakdur_std);
%         plot(xi,f,'Color',colmap(s.lig_dose == doses,:))
%         set(gca,'XLim',[-20 100])
%         title('Peak dur std')
%     end
    
    if ~isempty(peakdis_mean)
        subplot(nrows,ncols,ncols*(reorder(s.lig_index)-1)+7)
        hold on
        plot(find(s.lig_dose == doses),median(peakdis_mean),'x')
        set(gca,'XLim',[.5 length(doses)+.5],'YLim',[50 160])
        set(gca,'XTick',1:length(doses),'XTickLabel',doses)
        title('Peak distance')
    end
    
%     if sum(~isnan(peakdis_std))
%         subplot(nrows,ncols,ncols*(reorder(s.lig_index)-1)+9)
%         hold on
%         [f,xi] = ksdensity(peakdis_std);
%         plot(xi,f,'Color',colmap(s.lig_dose == doses,:))
%         set(gca,'XLim',[-20 100])
%         title('Peak dis std')
%     end
    
end

return

%%
close all

f4 = figure;

xfac = 1;
yfac = 1;
fontsize = 4;

setFigure(f4,xfac,yfac,fontsize)

nrows = 7;
ncols = 9; % nEdge; SNR; ...

for isite = sites_all
    
    [radial_dists c_signal_tmp tmp2 nEdge SNR amp pw peakdur_mean peakdur_std peakdis_mean peakdis_std] = edge_snr_score_pw_distdur(isite);
    
    s = siteprop(isite);
    
    if ~isempty(radial_dists)
        subplot(nrows,ncols,ncols*(s.lig_index-1)+1)
        hold on
        [f,xi] = ksdensity(radial_dists,'width',0.1);
        plot(xi,f,'Color',colmap(s.lig_dose == doses,:))
        set(gca,'XLim',[-.2 1.5])
        ylabel(s.lig_name)
        title('Final score')
    end
    
    if ~isempty(nEdge)
        subplot(nrows,ncols,ncols*(s.lig_index-1)+2)
        hold on
        [f,xi] = ksdensity(nEdge,'width',1);
        plot(xi,f,'Color',colmap(s.lig_dose == doses,:))
        set(gca,'XLim',[-2 10])
        title('nEdge')
    end
    
    if ~isempty(SNR)
        subplot(nrows,ncols,ncols*(s.lig_index-1)+3)
        hold on
        [f,xi] = ksdensity(SNR);
        plot(xi,f,'Color',colmap(s.lig_dose == doses,:))
        set(gca,'XLim',[-5 30])
        title('SNR')
    end
    
    if ~isempty(amp)
        subplot(nrows,ncols,ncols*(s.lig_index-1)+4)
        hold on
        [f,xi] = ksdensity(amp);
        plot(xi,f,'Color',colmap(s.lig_dose == doses,:))
        set(gca,'XLim',[-0.002 0.03])
        title('Amplitude')
    end
    
    if ~isempty(pw)
        subplot(nrows,ncols,ncols*(s.lig_index-1)+5)
        hold on
        [f,xi] = ksdensity(pw);
        plot(xi,f,'Color',colmap(s.lig_dose == doses,:))
        set(gca,'XLim',[-0.01 .1])
        title('Pairwise')
    end
    
    if ~isempty(peakdur_mean)
        subplot(nrows,ncols,ncols*(s.lig_index-1)+6)
        hold on
        [f,xi] = ksdensity(peakdur_mean);
        plot(xi,f,'Color',colmap(s.lig_dose == doses,:))
        set(gca,'XLim',[0 300])
        title('Peak duration')
    end
    
    if sum(~isnan(peakdur_std))
        subplot(nrows,ncols,ncols*(s.lig_index-1)+7)
        hold on
        [f,xi] = ksdensity(peakdur_std);
        plot(xi,f,'Color',colmap(s.lig_dose == doses,:))
        set(gca,'XLim',[-20 100])
        title('Peak dur std')
    end
    
    if ~isempty(peakdis_mean)
        subplot(nrows,ncols,ncols*(s.lig_index-1)+8)
        hold on
        [f,xi] = ksdensity(peakdis_mean);
        plot(xi,f,'Color',colmap(s.lig_dose == doses,:))
        set(gca,'XLim',[0 300])
        title('Peak distance')
    end
    
    if sum(~isnan(peakdis_std))
        subplot(nrows,ncols,ncols*(s.lig_index-1)+9)
        hold on
        [f,xi] = ksdensity(peakdis_std);
        plot(xi,f,'Color',colmap(s.lig_dose == doses,:))
        set(gca,'XLim',[-20 100])
        title('Peak dis std')
    end
    
end
