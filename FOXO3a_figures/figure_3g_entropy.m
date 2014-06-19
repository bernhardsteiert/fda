% Fig. 3G: Mutual information of early PCs and pulsing decision
load('./Workspaces/mi_boot_binary')

rng(0)
n = 10000;

legh = nan(1,3);
legstr = {'MI bootstrap','MI early PCs','Entropy Pulsing'};
f1 = figure;
hold on

xfac = 1;
yfac = .85;
fontsize = 16;
setFigure(f1,xfac,yfac,fontsize)

hist(mi_boot,15)
ylim = get(gca,'YLim');
mean_boot = mean(mi_boot);
std_boot = std(mi_boot);
legh(1) = plot([mean_boot-3*std_boot mean_boot+3*std_boot],.95*ylim(2)*[1 1],'LineWidth',2);
plot([mean_boot-3*std_boot mean_boot-3*std_boot],ylim,':')
plot([mean_boot+3*std_boot mean_boot+3*std_boot],ylim,':')
xlim = get(gca,'XLim');
legh(2) = plot(2*[xlim(2) xlim(2)],ylim,'k','LineWidth',2);
xlim2 = get(gca,'XLim');
legh(3) = plot(1.6*[xlim2(2) xlim2(2)],ylim,'r','LineWidth',2);
set(gca,'XTick',[0 mean_boot+3*std_boot 2*xlim(2) 1.6*xlim2(2)],'XTickLabel',[0 mean_boot+3*std_boot mi_pcs_dists h_dists])
xt = range(get(gca,'XLim'))/35;
yt = ylim(2)/15;
breakpos = 3*std_boot + (2*xlim(2)-3*std_boot)/2;
plot([breakpos-xt breakpos+xt]-xt/3,[0 yt],'Color',[.7 .7 .7])
plot([breakpos-xt breakpos+xt]+xt/3,[0 yt],'Color',[.7 .7 .7])
breakpos2 = 2*xlim(2) + (1.6*xlim2(2)-2*xlim(2))/2;
plot([breakpos2-xt breakpos2+xt]-xt/3,[0 yt],'Color',[.7 .7 .7])
plot([breakpos2-xt breakpos2+xt]+xt/3,[0 yt],'Color',[.7 .7 .7])
legend(legh,legstr)

xlabel('Information [bit]')
set(gca,'YTick',[])
xlim3 = get(gca,'XLim');
h = axes('Position',get(gca,'Position'));
set(h,'XAxisLocation','top','Color','None','YTick',[],'XLim',xlim3,'XTick',[0 mean_boot+3*std_boot 2*xlim(2) 1.6*xlim2(2)],'XTickLabel',{'0', sprintf('%3.3g',(mean_boot+3*std_boot)/h_dists*100), sprintf('%3.3g',mi_pcs_dists/h_dists*100), sprintf('%3.3g',h_dists/h_dists*100)})
xlabel('Information [%]')