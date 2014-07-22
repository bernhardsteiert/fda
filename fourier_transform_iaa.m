%% Analyse data from Jerry
c_signal_single = {};
c_signal_single{1} = csvread('./IAA/output.ratio-s1.csv')';
c_signal_single{2} = csvread('./IAA/output.ratio-s2.csv')';
c_signal_single{3} = csvread('./IAA/output.ratio-s3.csv')';
c_signal_single{4} = csvread('./IAA/output.ratio-s4.csv')';

timestamp = 0:5:5*size(c_signal_single{1},1)-1;
timestamp = timestamp';

% labels = {'< 10%', '10-25%', '25-50%', '50-75%', '75-90%', '> 90%', 'Pink noise', 'Sine + white noise'};
labels = {'output.ratio-s1','output.ratio-s2','output.ratio-s3','output.ratio-s4'};
colmap = hsv(length(c_signal_single));
% colmap(1,:) = [0 0 0];

figure

linewidth = 2;

xfac = 1;
yfac = .7;
fontsize = 10;

setFigure(gcf,xfac,yfac,fontsize)
hold on

range_ind = 1:length(timestamp);

Fs = 1./((timestamp(2)-timestamp(1))*60); % Sampling every 5 min
L = length(range_ind);

NFFT = 2^nextpow2(L);

legh = [];
legstr = {};

% for iplot = [1:6 8 7]
for iplot = 1:4

    c_signal_fft = [];
    c_signal_autocorr = [];
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
        
        c_signal_autocorr = [c_signal_autocorr autom(c_signal_tmp(range_ind,i))'];
    end
    
    c_signal_autocorr_norm = c_signal_autocorr./repmat(c_signal_autocorr(1,:),size(c_signal_autocorr,1),1);
    mean_auto = nanmean(c_signal_autocorr_norm,2);
    std_auto = nanstd(c_signal_autocorr_norm,[],2);
    std_auto = std_auto./sqrt(size(c_signal_autocorr_norm,2)-1); % Std of mean --> 1/sqrt(N-1)
    
    mean_fft = nanmean(c_signal_fft,2);
    std_fft = nanstd(c_signal_fft,[],2);
    std_fft = std_fft./sqrt(size(c_signal_fft,2)-1); % Std of mean --> 1/sqrt(N-1)
    
    subplot(1,2,1)
    hold on
    
    legh = [legh semilogx(f,10*log10(mean_fft),'Color',colmap(iplot,:),'LineWidth',linewidth)];
    
    tmpx = [f'; flipud(f')];
    % 3 sigma == 99.7% confidence 
    tmpy = [10*log10(mean_fft+3*std_fft./2); flipud(10*log10(mean_fft-3*std_fft./2))];

    ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp, 'FaceColor', colmap(iplot,:)*0.1+0.9, 'EdgeColor', colmap(iplot,:)*0.1+0.9,'LineWidth',linewidth);
    ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', colmap(iplot,:)*0.3+0.7,'LineWidth',linewidth);
    
    subplot(1,2,2)
    hold on
    
    plot(timestamp(1:size(c_signal_autocorr,1)),mean_auto,'Color',colmap(iplot,:),'LineWidth',linewidth)
    
    tmpx = [timestamp(1:size(c_signal_autocorr,1)); flipud(timestamp(1:size(c_signal_autocorr,1)))];
    % 3 sigma == 99.7% confidence 
    tmpy = [mean_auto+3*std_auto./2; flipud(mean_auto-3*std_auto./2)];

    ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp, 'FaceColor', colmap(iplot,:)*0.1+0.9, 'EdgeColor', colmap(iplot,:)*0.1+0.9,'LineWidth',linewidth);
    ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', colmap(iplot,:)*0.3+0.7,'LineWidth',linewidth);
    
%     legstr{end+1} = labels{iplot};

end

subplot(1,2,1)
set(gca,'XLim',[5e-5 1.2e-3])
set(gca,'XScale','log');
% set(gca,'YLim',[-95 -55])

% plot([1/(80*60) 1/(80*60)],get(gca,'YLim'),'k:') % equals 85min period
plot([1/(300*60) 1/(300*60)],get(gca,'YLim'),'k:','LineWidth',linewidth) % equals 300min period
plot([1/(15*60) 1/(15*60)],get(gca,'YLim'),'k:','LineWidth',linewidth) % equals 15min period

% set(gca,'XTick',[1/(300*60) 1/(80*60) 1/(15*60)],'XTickLabel',{'1/(300*60)','1/(80*60)','1/(15*60)'})
% title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('log10 |Y(f)|')
ylabel('|Y(f)|^2 [dB]')

legend(legh,labels,'Location','SouthWest')

xlim3 = get(gca,'XLim');
h = axes('Position',get(gca,'Position'));
set(h,'XAxisLocation','top','Color','None','YTick',[],'XLim',xlim3,'XTick',[1/(300*60) 1/(80*60) 1/(15*60)],'XTickLabel',{'300','80','15'},'XScale','log')
xlabel('Period [min]')

subplot(1,2,2)

title('Autocorrelation')
xlabel('time [min]')
ylabel('ACF')
set(gca,'XLim',[0 50])
