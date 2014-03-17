%% Fourier Analysis of c_signal
% http://www.wavemetrics.com/products/igorpro/dataanalysis/signalprocessing/powerspectra.htm:
% Power spectra can be computed for the entire signal at once (a "periodogram") or periodograms of segments of the time signal can be averaged together to form the "power spectral density".

% % Generate workspace
% 
% close all
% clear all
% clc
% 
% addpath('./noise_generation/')
% 
% remotepath = mypath();
% 
% fdaMPath = [remotepath 'fda'];
% addpath(fdaMPath)
% 
% grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
% addpath(grabdataPath)
% 
% plot_sites = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:70];
% lignames = {'EGF','IGF','FGF','HRG','HGF','EPR','BTC'};
% 
% times = cell(0);
% signals = cell(0);
% signals_raw = cell(0);
% celltype = [];
% 
% for isite = plot_sites
%     if exist(remotepath,'dir')
%         [times{end+1},intensity] = grabdata(isite);
%     else
%         load(['./Workspaces/site_' num2str(isite)])
%         times{end+1} = timestamp;
%     end
% 
%     log_trafo = 1; % log-transform signal
% 
%     if log_trafo
%         signals_raw{end+1} = log10(intensity);
%     else
%         signals_raw{end+1} = intensity;
%     end
%     
%     signals{end+1} = signals_raw{end} - repmat(nanmean(signals_raw{end},2),1,size(signals_raw{end},2));
%     
%     celltype = [celltype ones(1,size(intensity,2))*isite];
% end
% 
% timestamp = times{1}; % same time sampling for all data sets
% c_signal = cell2mat(signals);
% c_signal_raw = cell2mat(signals_raw);
% 
% dists = [];
% celltypeharm = [];
% 
% clear radial_dist
% 
% allligs = zeros(size(celltype));
% 
% for isite = plot_sites
%     sprop = siteprop(isite);
%     allligs = allligs + sprop.lig_index*(celltype==isite);
%     
%     [radial_dists c_signal_tmp tmp2 nEdge SNR amp pw peakdur_mean peakdur_std peakdis_mean peakdis_std] = edge_snr_score_pw_distdur(isite);
%     
%     dists = [dists radial_dists];
%     
%     celltypeharm = [celltypeharm ones(size(radial_dists))*isite];
%     
% end
% 
% [radial_dist_sorted ind_sort_radial] = sort(dists);
% 
% 
% nsignal = 40; % #signals to be used for mean calculation
% % for ip = 1:length(plot_sites)
% thres_sorted = .1;
% ind_sort_radial_top = ind_sort_radial(end-round(length(ind_sort_radial)*thres_sorted)+1:end);
% ind_sort_radial_bottom = ind_sort_radial(1:round(length(ind_sort_radial)*thres_sorted));
% 
% c_signal_single = {};
% c_signal_single{1} = c_signal;
% c_signal_single{2} = c_signal(:,ind_sort_radial_top);
% c_signal_single{3} = c_signal(:,ind_sort_radial_bottom);
% 
% % White noise:
% c_signal_single{4} = .002*repmat(sin(2*pi*timestamp*60/(80*60)),1,size(c_signal,2))+.001*randn(size(c_signal));
% % Pink noise:
% c_signal_single{5} = c_signal;
% c_signal_single{5}(:) = .01*pinknoise(length(c_signal(:)))';
% 
% save('fourier_signals','timestamp','c_signal_single','thres_sorted')

% return

%% Plot
load('./Workspaces/fourier_signals')

labels = {'Main dataset', sprintf('Only %i%% highest scores',thres_sorted*100), sprintf('Only %i%% lowest scores',thres_sorted*100), 'White noise \\w Sinusoidal', 'Pink noise'};
colmap = lines(length(labels));
colmap(1,:) = [0 0 0];

figure

linewidth = 1;

hold on

time_range = [200 500];
[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

Fs = 1./((timestamp(2)-timestamp(1))*60); % Sampling every 5 min
L = length(range_ind);

NFFT = 2^nextpow2(L);

legh = [];
legstr = {};

for iplot = 1:length(c_signal_single)

    c_signal_fft = [];
    for i = 1:size(c_signal_single{iplot},2)
        Y = fft(c_signal_single{iplot}(range_ind,end-i+1),NFFT)/L;

        f = Fs/2*linspace(0,1,NFFT/2+1);

        c_signal_fft = [c_signal_fft 2*abs(Y(1:NFFT/2+1)).^2];
    end
    
    mean_fft = mean(c_signal_fft,2);
    std_fft = std(c_signal_fft,[],2);
    std_fft = std_fft./sqrt(size(c_signal_fft,2)-1); % Std of mean --> 1/sqrt(N-1)
    
    legh = [legh semilogx(f,10*log10(mean_fft),'Color',colmap(iplot,:),'LineWidth',linewidth)];
    
    tmpx = [f'; flipud(f')];
    % 3 sigma == 99.7% confidence 
    tmpy = [10*log10(mean_fft+3*std_fft./2); flipud(10*log10(mean_fft-3*std_fft./2))];

    ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp, 'FaceColor', colmap(iplot,:)*0.1+0.9, 'EdgeColor', colmap(iplot,:)*0.1+0.9,'LineWidth',linewidth);
    ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', colmap(iplot,:)*0.3+0.7,'LineWidth',linewidth);
    
    legstr{end+1} = labels{iplot};

end


set(gca,'XLim',[5e-5 1.2e-3])
set(gca,'XScale','log');
set(gca,'YLim',[-80 -50])

% plot([1/(80*60) 1/(80*60)],get(gca,'YLim'),'k:') % equals 85min period
plot([1/(300*60) 1/(300*60)],get(gca,'YLim'),'k:','LineWidth',linewidth) % equals 300min period
plot([1/(15*60) 1/(15*60)],get(gca,'YLim'),'k:','LineWidth',linewidth) % equals 15min period

% set(gca,'XTick',[1/(300*60) 1/(80*60) 1/(15*60)],'XTickLabel',{'1/(300*60)','1/(80*60)','1/(15*60)'})
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('log10 |Y(f)|')
ylabel('|Y(f)|^2 [dB]')

legend(legh,legstr)