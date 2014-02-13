%% Fourier Analysis of c_signal

close all
% clear all
clc

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
addpath(grabdataPath)


% plot_sites = [10 17 64];
% plot_name = {'No Lig','IGF','BTC'};
% plot_sites = [64];
% plot_sites = [37 64]; % HRG and BTC
plot_name = {'HRG','BTC'};

% plot_sites = [65 64]; % HRG and BTC
% plot_name = {'BTC 50','BTC 100'};

% plot_sites = [17 64]; % IGF and BTC
% plot_name = {'HRG','BTC'};

plot_sites = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:70];

times = cell(0);
signals = cell(0);
signals_raw = cell(0);
celltype = [];

for isite = plot_sites
    if exist(remotepath,'dir')
        [times{end+1},intensity] = grabdata(isite);
    else
        load(['./Workspaces/site_' num2str(isite)])
        times{end+1} = timestamp;
    end

    log_trafo = 1; % log-transform signal

    if log_trafo
        signals_raw{end+1} = log10(intensity);
    else
        signals_raw{end+1} = intensity;
    end
    
    signals{end+1} = signals_raw{end} - repmat(nanmean(signals_raw{end},2),1,size(signals_raw{end},2));
    
    celltype = [celltype ones(1,size(intensity,2))*isite];
end

timestamp = times{1}; % same time sampling for all data sets
c_signal = cell2mat(signals);
c_signal_raw = cell2mat(signals_raw);

dists = [];
celltypeharm = [];

nEdges = [];
SNRs = [];
amps = [];
pws = [];
c_signal_woNharm = [];

peakdur_means = [];
peakdur_stds = [];
peakdis_means = [];
peakdis_stds = [];

clear radial_dist

for isite = plot_sites
    
%     radial_dists = radial_dist(isite);
%     radial_dists = freq_analysis(isite);
%     radial_dists = edge_snr_score(isite);
    [radial_dists c_signal_tmp tmp2 nEdge SNR amp pw peakdur_mean peakdur_std peakdis_mean peakdis_std] = edge_snr_score_pw_distdur(isite);
%     radial_dists = radial_dist_pw(isite);
    
    dists = [dists radial_dists];
    
    celltypeharm = [celltypeharm ones(size(radial_dists))*isite];
    
    nEdges = [nEdges nEdge];
    SNRs = [SNRs SNR];
    amps = [amps amp];
    pws = [pws pw];
%     c_signal_woNharm = [c_signal_woNharm c_signal_tmp];
    
    peakdur_means = [peakdur_means peakdur_mean];
    peakdur_stds = [peakdur_stds peakdur_std];
    peakdis_means = [peakdis_means peakdis_mean];
    peakdis_stds = [peakdis_stds peakdis_std];
    
end

[radial_dist_sorted ind_sort_radial] = sort(dists);
% peakdis_means(peakdis_means>290) = 0;
% [radial_dist_sorted ind_sort_radial] = sort(peakdis_means); % For checking IGF

%%
f = figure;

xfac = 2;
yfac = 1;
fontsize = 12;

setFigure(f,xfac,yfac,fontsize)

dist_thres = .1;
valid_dist_inds = find(dists > dist_thres);

tested_score = nEdges; % SNRs or nEdges or anything else
baredges = linspace(min(tested_score),max(tested_score),9);

% tested_score = SNRs; % SNRs or nEdges or anything else
baredges = linspace(min(tested_score),max(tested_score),9);

% tested_score = amps; % SNRs or nEdges or anything else
% tested_score = pws; % SNRs or nEdges or anything else
% baredges = linspace(min(tested_score),max(tested_score),10);

nrows = length(plot_sites);
ncols = 9;

for ip = 1:length(plot_sites)
    tested_inds = find(dists > dist_thres & celltypeharm == plot_sites(ip));
    
    subplot(nrows,ncols,(ip-1)*ncols+1)
    
    baredges = linspace(min(dists),max(dists),9);
    bar(baredges,histc(dists(tested_inds),baredges));
    title('PulsStr')
    ylabel(plot_name{ip})
    set(gca,'XLim',[min(dists)-0.02 max(dists)]*1.1)
    
    subplot(nrows,ncols,(ip-1)*ncols+2)
    
    baredges = linspace(min(nEdges),max(nEdges),9);
    bar(baredges,histc(nEdges(tested_inds),baredges));
    title('nEdge')
    set(gca,'XLim',[min(nEdges)+1 max(nEdges)]*1.1)
    
    subplot(nrows,ncols,(ip-1)*ncols+3)
    
    baredges = linspace(min(SNRs),max(SNRs),9);
    bar(baredges,histc(SNRs(tested_inds),baredges));
    title('SNR')
    set(gca,'XLim',[min(SNRs)-5 max(SNRs)]*1.1)
    
    subplot(nrows,ncols,(ip-1)*ncols+4)
    
    baredges = linspace(min(amps),max(amps),9);
    bar(baredges,histc(amps(tested_inds),baredges));
    title('Amplitude')
    set(gca,'XLim',[min(amps)-0.005 max(amps)]*1.1)
    
    subplot(nrows,ncols,(ip-1)*ncols+5)
    
    baredges = linspace(min(pws),max(pws),9);
    bar(baredges,histc(pws(tested_inds),baredges));
    title('Pairwise')
    set(gca,'XLim',[min(pws) max(pws)]*1.1)
    
    subplot(nrows,ncols,(ip-1)*ncols+6)
    
    baredges = linspace(min(peakdur_means),150,9);
    bar(baredges,histc(peakdur_means(tested_inds),baredges));
    title('Peak duration (min)')
    set(gca,'XLim',[min(peakdur_means)-20 150]*1.1)
    
    subplot(nrows,ncols,(ip-1)*ncols+7)
    
    baredges = linspace(min(peakdur_stds),max(peakdur_stds),9);
    bar(baredges,histc(peakdur_stds(tested_inds),baredges));
    title('Peak dur std')
    set(gca,'XLim',[min(peakdur_stds)-10 max(peakdur_stds)]*1.1)
    
    subplot(nrows,ncols,(ip-1)*ncols+8)
    
    baredges = linspace(min(peakdis_means),150,9);
    bar(baredges,histc(peakdis_means(tested_inds),baredges));
    title('Peak distance')
    set(gca,'XLim',[min(peakdis_means)-20 150]*1.1)
    
    subplot(nrows,ncols,(ip-1)*ncols+9)
    
    baredges = linspace(min(peakdis_stds),max(peakdis_stds),9);
    bar(baredges,histc(peakdis_stds(tested_inds),baredges));
    title('Peak dis std')
    set(gca,'XLim',[min(peakdis_stds)-10 max(peakdis_stds)]*1.1)
    
end

% return

for ip = 1:length(plot_sites)
    isite = plot_sites(ip);
    figure
    
    posFig = get(gcf,'Position');
%     posFig(3) = posFig(3)/2;
%     posFig(4) = posFig(4)/2;
    set(gcf,'Position',posFig)
    set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./15);
    
%     c_signal_single = c_signal(:,ind_sort_radial);
    c_signal_single = c_signal_raw(:,ind_sort_radial);
%     c_signal_single = c_signal_woNharm(:,ind_sort_radial);
    c_signal_single = c_signal_single(:,(celltypeharm(ind_sort_radial) == plot_sites(ip)));
    % Remove time mean
%     c_signal_single = c_signal_single - repmat(mean(c_signal_single,1),size(c_signal_single,1),1);
    
    time_range = [200 500];
    [tmp range_ind_min] = min(abs(timestamp - time_range(1)));
    [tmp range_ind_max] = min(abs(timestamp - time_range(2)));
    range_ind = range_ind_min:range_ind_max;
    
    nrows = 5;
    ncols = 4;
    
    mydists = radial_dist_sorted(celltypeharm(ind_sort_radial) == plot_sites(ip));
    
    for isig = 1:(nrows*ncols)
        subplot(nrows,ncols,isig)
        
        plot(timestamp(range_ind),c_signal_single(range_ind,end-isig+1)) % For c_signal and c_signal_raw
%         plot(timestamp(range_ind),c_signal_single(range_ind,21-isig))

%         plot(timestamp(tmp2),c_signal_single(:,end-isig+1))
        
        set(gca,'XLim',time_range)
        set(gca,'YLim',[-.02 .03])
%         set(gca,'YLim',[.95 1.1])

        title(sprintf('PulsStr = %g',mydists(end-isig+1)))
        
        if isig == (nrows-1)*ncols+1
            xlabel('time [min]')
        end
        if isig == 1
            ylabel('log_{10} FOXO3a Cyt/Nuc ratio');
        end
    end
    
end

return

%%
for ip = 1:length(plot_sites)
    isite = plot_sites(ip);
    figure
    
    posFig = get(gcf,'Position');
%     posFig(3) = posFig(3)/2;
%     posFig(4) = posFig(4)/2;
    set(gcf,'Position',posFig)
    set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./15);
    
    c_signal_single = c_signal(:,ind_sort_radial);
    c_signal_single = c_signal_single(:,(celltypeharm(ind_sort_radial) == plot_sites(ip)));
    % Remove time mean
    c_signal_single = c_signal_single - repmat(mean(c_signal_single,1),size(c_signal_single,1),1);

    time_range = [200 500];
    [tmp range_ind_min] = min(abs(timestamp - time_range(1)));
    [tmp range_ind_max] = min(abs(timestamp - time_range(2)));
    range_ind = range_ind_min:range_ind_max;

    nrows = 5;
    ncols = 4;
    
    Fs = 1./((timestamp(2)-timestamp(1))*60); % Sampling every 5 min
    L = length(range_ind);

    NFFT = 2^nextpow2(L);
    for i = 1:(nrows*ncols)
        subplot(nrows,ncols,i)
        
        Y = fft(c_signal_single(range_ind,end-i+1),NFFT)/L;

        f = Fs/2*linspace(0,1,NFFT/2+1);
        plot(f,2*abs(Y(1:NFFT/2+1)))
        
        hold on
        
        set(gca,'YLim',[0 0.012])
        plot([2e-4 2e-4],get(gca,'YLim'),'k:') % equals 83min period
       
        if i == 1
            title('Single-Sided Amplitude Spectrum of y(t)')
            ylabel('|Y(f)|')
        end
        if i == (nrows-1)*ncols+1
            xlabel('Frequency (Hz)')
        end
    end
end


%%

nsignal = 40; % #signals to be used for mean calculation
% for ip = 1:length(plot_sites)
ind_sort_radial_top1000 = ind_sort_radial(4280:end);
ip = 1;
    isite = plot_sites(ip);
    figure
    
    c_signal_single = c_signal(:,ind_sort_radial_top1000);
%     c_signal_single = c_signal_single(:,(celltypeharm(ind_sort_radial) == plot_sites(ip)));

    time_range = [200 500];
    [tmp range_ind_min] = min(abs(timestamp - time_range(1)));
    [tmp range_ind_max] = min(abs(timestamp - time_range(2)));
    range_ind = range_ind_min:range_ind_max;
    
%     c_signal_single = c_signal_single-repmat(mean(c_signal_single(range_ind,:),1),size(c_signal_single,1),1);
    
    Fs = 1./((timestamp(2)-timestamp(1))*60); % Sampling every 5 min
    L = length(range_ind);

    NFFT = 2^nextpow2(L);
    
    c_signal_fft = [];
%     for i = 1:nsignal
    for i = 1:size(c_signal_single,2)
        Y = fft(c_signal_single(range_ind,end-i+1),NFFT)/L;

        f = Fs/2*linspace(0,1,NFFT/2+1);
        
%         if Y(5) > mean([Y(4) Y(6)])
            c_signal_fft = [c_signal_fft 2*abs(Y(1:NFFT/2+1))];
%         end
    end
    
    plot(f,log10(mean(c_signal_fft,2)))
    hold on
    plot([2e-4 2e-4],get(gca,'YLim'),'k:') % equals 83min period
%     [fi xi] = ksdensity(mean(c_signal_fft,2));
%     plot(xi,fi)
%     title(sprintf('Single-Sided Amplitude Spectrum of y(t); Mean spectrum of %i signals',nsignal))
    title(sprintf('Single-Sided Amplitude Spectrum of y(t); Mean spectrum of %i signals',size(c_signal_single,2)))
    xlabel('Frequency (Hz)')
    ylabel('log10 |Y(f)|')
% end