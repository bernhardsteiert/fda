%% Fourier Analysis of c_signal
% http://www.wavemetrics.com/products/igorpro/dataanalysis/signalprocessing/powerspectra.htm:
% Power spectra can be computed for the entire signal at once (a "periodogram") or periodograms of segments of the time signal can be averaged together to form the "power spectral density".

% Generate workspace

% close all
% clear all
% clc
% 
% myextension = '130722_corrected_retracked_all_cleaned';
% timeshift = 0;
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
% plot_sites = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:69];
% lignames = {'EGF','IGF','FGF','HRG','HGF','EPR','BTC'};
% 
% times = cell(0);
% signals = cell(0);
% signals_raw = cell(0);
% celltype = [];
% 
% for isite = plot_sites
% %     if exist(remotepath,'dir')
% %         [times{end+1},intensity] = grabdata(isite);
% %     else
%         load(['./Workspaces/site_' num2str(isite) '_' myextension])
%         times{end+1} = timestamp;
% %     end
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
% 
% for isite = plot_sites
%     [radial_dists tmp1 tmp2 nEdge SNR] = edge_snr_score_pw_distdur(isite,myextension,timeshift,1/120,'harm_basis_130722_corrected_retracked_all_cleaned_late',1);
%     dists = [dists radial_dists];
% end
% 
% [radial_dist_sorted ind_sort_radial] = sort(dists);
% 
% thres_sorted = [0 .1 .25 .5 .75 .9 1];
% percentile_ind = thres_sorted * round(length(dists));
% 
% ind_sort_radial_top = ind_sort_radial(end-round(length(ind_sort_radial)*thres_sorted)+1:end);
% ind_sort_radial_bottom = ind_sort_radial(1:round(length(ind_sort_radial)*thres_sorted));
% 
% c_signal_single = {};
% for ithres = 1:length(thres_sorted)-1
%     myind = ind_sort_radial(percentile_ind(ithres)+1:percentile_ind(ithres+1));
%     c_signal_single{end+1} = c_signal(:,myind);
% end
% 
% % White noise:
% c_signal_single{end+1} = .002*repmat(sin(2*pi*timestamp*60/(80*60)),1,size(c_signal,2))+.001*randn(size(c_signal));
% % Pink noise:
% c_signal_single{end+1} = c_signal;
% c_signal_single{end}(:) = .01*pinknoise(length(c_signal(:)))';
% % c_signal_single{6} = c_signal;
% % c_signal_single{6}(:) = .01*pinknoise2(length(c_signal(:)))';
% 
% save('fourier_signals_corrected_cleaned_newBTC','timestamp','c_signal_single','celltype')
% 
% return

%% Plot
addpath('./Functions/')

thres_sorted = [0 .1 .25 .5 .75 .9 1];
% load('fourier_signals')
% time_range = [200 500];
load('./Workspaces/fourier_signals_corrected_cleaned_newBTC')
time_range = [200 1475];

labels = {'< 10%', '10-25%', '25-50%', '50-75%', '75-90%', '> 90%', 'Pink noise', 'Sine + white noise'};
colmap = hsv(length(c_signal_single));
% colmap(1,:) = [0 0 0];

figure

linewidth = 2;

xfac = 1;
yfac = 1;
fontsize = 16;

setFigure(gcf,xfac,yfac,fontsize)
hold on

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

Fs = 1./((timestamp(2)-timestamp(1))*60); % Sampling every 5 min
L = length(range_ind);

NFFT = 2^nextpow2(L);

legh = [];
legstr = {};

c_signal_single{8} = c_signal_single{8}*.6;

for iplot = [1:6 8 7]

    c_signal_fft = [];
    c_signal_tmp = c_signal_single{iplot};
    c_signal_tmp(isinf(c_signal_tmp)) = nan;
    for i = 1:size(c_signal_tmp,2)
        if sum(isnan(c_signal_tmp(:,i))) > length(c_signal_tmp(:,i))-2
            c_signal_tmp(:,i) = 0;
        end
        if sum(isnan(c_signal_tmp(:,i))) > 0
            c_signal_tmp(:,i) = interp1(timestamp(~isnan(c_signal_tmp(:,i))),c_signal_tmp(~isnan(c_signal_tmp(:,i)),i),timestamp);
            vec = ~isnan(c_signal_tmp(:,i))';
            rl = find(vec ~= [vec(2:end), vec(end)+1]);
            data =  vec(rl);
            rl(2:end) = rl(2:end) - rl(1:end-1);
            if ~data(1)
                c_signal_tmp(1:rl(1),i) = c_signal_tmp(rl(1)+1,i);
            end
            if ~data(end)
                c_signal_tmp(end-rl(end)+1:end,i) = c_signal_tmp(end-rl(end),i);
            end
        end
            
        Y = fft(triang(length(range_ind)).*c_signal_tmp(range_ind,i),NFFT)/L;

        f = Fs/2*linspace(0,1,NFFT/2+1);

        c_signal_fft = [c_signal_fft 2*abs(Y(1:NFFT/2+1)).^2];
    end
    
    mean_fft = nanmean(c_signal_fft,2);
    std_fft = nanstd(c_signal_fft,[],2);
    std_fft = std_fft./sqrt(size(c_signal_fft,2)-1); % Std of mean --> 1/sqrt(N-1)
    
    legh = [legh semilogx(f,10*log10(mean_fft),'Color',colmap(iplot,:),'LineWidth',linewidth)];
    
    tmpx = [f'; flipud(f')];
    % 3 sigma == 99.7% confidence 
    tmpy = [10*log10(mean_fft+3*std_fft./2); flipud(10*log10(mean_fft-3*std_fft./2))];

    ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp, 'FaceColor', colmap(iplot,:)*0.1+0.9, 'EdgeColor', colmap(iplot,:)*0.1+0.9,'LineWidth',linewidth);
    ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', colmap(iplot,:)*0.3+0.7,'LineWidth',linewidth);
    
%     legstr{end+1} = labels{iplot};

end


set(gca,'XLim',[5e-5 1.2e-3])
set(gca,'XScale','log');
set(gca,'YLim',[-95 -60])

% plot([1/(80*60) 1/(80*60)],get(gca,'YLim'),'k:') % equals 85min period
plot([1/(300*60) 1/(300*60)],get(gca,'YLim'),'k:','LineWidth',linewidth) % equals 300min period
plot([1/(15*60) 1/(15*60)],get(gca,'YLim'),'k:','LineWidth',linewidth) % equals 15min period

% set(gca,'XTick',[1/(300*60) 1/(80*60) 1/(15*60)],'XTickLabel',{'1/(300*60)','1/(80*60)','1/(15*60)'})
% title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('log10 |Y(f)|')
ylabel('|Y(f)|^2 [dB]')

legend(legh,labels)

xlim3 = get(gca,'XLim');
h = axes('Position',get(gca,'Position'));
set(h,'XAxisLocation','top','Color','None','YTick',[],'XLim',xlim3,'XTick',[1/(300*60) 1/(80*60) 1/(15*60)],'XTickLabel',{'300','80','15'},'XScale','log')
xlabel('Period [min]')
