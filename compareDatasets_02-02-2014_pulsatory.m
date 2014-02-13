extension = '02-02-2014';

sites_all = [49 48 25 24 1 50 47 26 23 2 51 46 27 22 3 45 28 21 4];
names_all = {'EGF 100','EGF 50','EGF 20','EGF 10'};

% sites_all = [53 44 29 20 5 54 43 30 19 6 55 42 31 18 7 56 41 32 17 8];
% names_all = {'BTC 100','BTC 50','BTC 20','BTC 10'};

% sites_all = [57 40 33 16 9 58 39 34 15 10 59 38 35 14 11 60 13 12];
% names_all = {'HRG 100','HRG 50','HRG 20','HRG 10'};

puls_thres = .5;

reorder = [1 2 3 4];

doses = [0 2.5/16 2.5/4 2.5 10];
colmap = flipud(winter(length(doses)));
ligdoses = [100 50 20 10];

nrows = 4;
% ncols = 9; % nEdge; SNR; ...
ncols = 8;

medians = nan(length(doses),nrows,ncols);
singles = [];
features = {'Final score','nEdge','SNR','Amplitude','Pairwise','Peak duration','Peak distance','Ratio pulsing cells'};
single_features = {'Final score','nEdge','SNR','Amplitude','Pairwise','Peak duration','Peak distance','Ratio pulsing cells'};

for isite = sites_all
    
    [radial_dists c_signal_tmp tmp2 nEdge SNR amp pw peakdur_mean peakdur_std peakdis_mean peakdis_std] = edge_snr_score_pw_distdur(isite,extension,0,1/180);
    
    s = siteprop(isite,extension);
    
    doseind = s.drug_dose == doses;
    medians(doseind,reorder(s.lig_dose == ligdoses),1) = median(radial_dists(radial_dists > puls_thres));
    medians(doseind,reorder(s.lig_dose == ligdoses),2) = median(nEdge(radial_dists > puls_thres));
    medians(doseind,reorder(s.lig_dose == ligdoses),3) = median(SNR(radial_dists > puls_thres));
    medians(doseind,reorder(s.lig_dose == ligdoses),4) = median(amp(radial_dists > puls_thres));
    medians(doseind,reorder(s.lig_dose == ligdoses),5) = median(pw(radial_dists > puls_thres));
    medians(doseind,reorder(s.lig_dose == ligdoses),6) = median(peakdur_mean(~isnan(peakdur_mean) & (radial_dists > puls_thres)));
    medians(doseind,reorder(s.lig_dose == ligdoses),7) = median(peakdis_mean(~isnan(peakdis_mean) & (radial_dists > puls_thres)));
    medians(doseind,reorder(s.lig_dose == ligdoses),8) = sum(radial_dists > puls_thres)/length(radial_dists);
    
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
            xlabel('MEKi dose')
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