% Monte-Carlo simulation of puls distance / puls width
close all

Nsample = 6137; % (size(singles,1))

addpath('./noise_generation/')

load('./Workspaces/site_1')
c_signal = log10(intensity);

time_range = [200 510];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

scale = mean(range(c_signal(range_ind,:),2));

wnoise = randn(length(range_ind),Nsample);
wnoise = wnoise*scale./4;
wnoise = wnoise.*repmat(rand(1,size(wnoise,2)),length(range_ind),1);

pnoise = nan(length(range_ind),Nsample);
for isample = 1:Nsample
    pnoise(:,isample) = pinknoise(length(range_ind));
end
pnoise = pnoise*scale./2;
pnoise = pnoise.*repmat(rand(1,size(pnoise,2)),length(range_ind),1);

subplot(1,2,1)
plot(timestamp(range_ind),wnoise(:,1))
set(gca,'XLim',time_range)
subplot(1,2,2)
plot(timestamp(range_ind),pnoise(:,1))
set(gca,'XLim',time_range)

[radial_dist c_signal_woNharm range_ind_white nEdges SNR amp pw peakdur_mean peakdur_std peakdis_mean peakdis_std] = edge_snr_score_pw_distdur_sim(timestamp(range_ind),wnoise);
[radial_dist_pink c_signal_woNharm_pink range_ind_pink nEdges_pink SNR_pink amp_pink pw_pink peakdur_mean_pink peakdur_std_pink peakdis_mean_pink peakdis_std_pink] = edge_snr_score_pw_distdur_sim(timestamp(range_ind),pnoise);

%%
close all
f1 = figure;

xfac = 1.5;
yfac = 1;
fontsize = 12;

setFigure(f1,xfac,yfac,fontsize)

xlim = [0 250];
ylim = [0 2.5];

baredges1 = linspace(0,250,26);
baredges2 = linspace(1.5,1.8,26);

subplot(3,3,1)
plot(peakdur_mean,radial_dist,'k.')
title('White Noise')
ylabel('Final score')
xlabel('Peak duration')
set(gca,'XLim',xlim,'YLim',ylim)
subplot(3,3,2)
plot(peakdur_mean_pink,radial_dist_pink,'k.')
title('Pink Noise')
set(gca,'XLim',xlim,'YLim',ylim)
subplot(3,3,3)
plot(singles(:,6),singles(:,1),'k.')
set(gca,'XLim',xlim,'YLim',ylim)
title('Main Dataset')
subplot(3,3,4)
bar(baredges1,histc(peakdur_mean,baredges1));
set(gca,'XLim',xlim)
subplot(3,3,5)
bar(baredges1,histc(peakdur_mean_pink,baredges1));
set(gca,'XLim',xlim)
title('Histograms of Peak duration')
subplot(3,3,6)
bar(baredges1,histc(singles(:,6),baredges1));
set(gca,'XLim',xlim)

[xp,lambda,c,Lmax] = boxcox(peakdur_mean_pink);
subplot(3,3,7)
bar(baredges2,histc(boxcox_apply(peakdur_mean,lambda,c),baredges2));
set(gca,'XLim',[min(baredges2) max(baredges2)])
subplot(3,3,8)
bar(baredges2,histc(xp,baredges2));
set(gca,'XLim',[min(baredges2) max(baredges2)])
title('After Box-Cox Trafo of Pink Noise distribution')
subplot(3,3,9)
bar(baredges2,histc(boxcox_apply(singles(:,6),lambda,c),baredges2));
set(gca,'XLim',[min(baredges2) max(baredges2)])