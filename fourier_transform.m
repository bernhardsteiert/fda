%% Fourier Analysis of c_signal

close all
clear all
clc

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
addpath(grabdataPath)


% plot_sites = [10 17 64];
% plot_name = {'No Lig','IGF','BTC'};
plot_sites = [64];
plot_name = {'BTC'};

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

clear radial_dist

for isite = plot_sites
    
    radial_dists = radial_dist(isite);
    
    dists = [dists radial_dists];
    
    celltypeharm = [celltypeharm ones(size(radial_dists))*isite];
    
end

[radial_dist_sorted ind_sort_radial] = sort(dists);


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
    c_signal_single = c_signal_single(:,(celltypeharm(ind_sort_radial) == plot_sites(ip)));
    % Remove time mean
    c_signal_single = c_signal_single - repmat(mean(c_signal_single,1),size(c_signal_single,1),1);
    
    time_range = [200 500];
    [tmp range_ind_min] = min(abs(timestamp - time_range(1)));
    [tmp range_ind_max] = min(abs(timestamp - time_range(2)));
    range_ind = range_ind_min:range_ind_max;
    
    nrows = 5;
    ncols = 4;
    
    for isig = 1:(nrows*ncols)
        subplot(nrows,ncols,isig)
        
        plot(timestamp(range_ind),c_signal_single(range_ind,end-isig+1))
%         plot(timestamp(range_ind),c_signal_single(range_ind,21-isig))
        
        set(gca,'XLim',time_range)
        set(gca,'YLim',[-.02 .02])
%         set(gca,'YLim',[.95 1.1])
        
        if isig == (nrows-1)*ncols+1
            xlabel('time [min]')
        end
        if isig == 1
            ylabel('log_{10} FOXO3a Cyt/Nuc ratio');
        end
    end
    
end


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

nsignal = 40; % #signals to be used for mean calculation
for ip = 1:length(plot_sites)
    isite = plot_sites(ip);
    figure
    
    c_signal_single = c_signal(:,ind_sort_radial);
    c_signal_single = c_signal_single(:,(celltypeharm(ind_sort_radial) == plot_sites(ip)));

    time_range = [200 500];
    [tmp range_ind_min] = min(abs(timestamp - time_range(1)));
    [tmp range_ind_max] = min(abs(timestamp - time_range(2)));
    range_ind = range_ind_min:range_ind_max;
    
    Fs = 1./((timestamp(2)-timestamp(1))*60); % Sampling every 5 min
    L = length(range_ind);

    NFFT = 2^nextpow2(L);
    
    c_signal_fft = [];
    for i = 1:nsignal
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
    title(sprintf('Single-Sided Amplitude Spectrum of y(t); Mean spectrum of %i signals',nsignal))
    xlabel('Frequency (Hz)')
    ylabel('log10 |Y(f)|')
end