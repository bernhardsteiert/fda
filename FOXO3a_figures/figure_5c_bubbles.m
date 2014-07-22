% Figure 5c: Bubble plot pulsing vs. mean response
addpath('./Functions/')
load('./Workspaces/bubble_04042014')

sites_all = [1:17 19:47 50:80];
extension = '04-04-2014_all_cleaned';

puls_thres = .4;

figure

subplot(1,3,1)
title('MCF10A')
hold on
subplot(1,3,2)
title('184A1')
hold on

marker = 'o';
legh = nan(2,2);

for isite = sites_all
    
    s = siteprop(isite,extension);
    if strmatch(s.celltype,'MCF10A','exact')
        isp = 1;
    else
        isp = 2;
    end
    subplot(1,3,isp)
    isp2 = 0;
    switch s.drug2_dose
        case 0
            markersize = 16;
        case .005
            markersize = 12;
            isp2 = 1;
        case .05
            markersize = 8;
        case .5
            markersize = 4;
    end
        
    if strmatch(s.lig_name,'IGF','exact')
        col1 = 1/3; % green
        ilig = 1;
    else
        col1 = 0; % red
        ilig = 2;
    end
    
    isp3 = 0;
    switch s.drug1_dose
        case 0
            col3 = .2;
        case 1/6^3
            col3 = .4;
        case 1/6^2
            col3 = .6;
        case 1/6
            col3 = .8;
        case 1
            col3 = 1;
            isp3 = 1;
    end
    
    color = hsv2rgb([col1 1 col3]);
    if isp2 && isp3
        legh(isp,ilig) = plot(mean(mean_response(celltype == isite)),sum(dists(celltype == isite) > puls_thres) / sum(celltype == isite),marker,'MarkerEdgeColor','k','MarkerSize',markersize,'MarkerFaceColor',color);
    else
        plot(mean(mean_response(celltype == isite)),sum(dists(celltype == isite) > puls_thres) / sum(celltype == isite),marker,'MarkerEdgeColor','k','MarkerSize',markersize,'MarkerFaceColor',color)
    end
    
    
end

subplot(1,3,1)
set(gca,'YLim',[-.1 .8])
xlabel('mean response (late PC1 - early PC1)')
ylabel('ratio pulsing cells')
legend(legh(1,:),{'IGF','EGF'},'Location','NorthWest')
subplot(1,3,2)
set(gca,'YLim',[-.1 .8])
xlabel('mean response (late PC1 - early PC1)')
ylabel('ratio pulsing cells')
legend(legh(2,:),{'IGF','EGF'},'Location','NorthWest')
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');

subplot(1,3,3)
hold on
set(gca,'XLim',xlim)
set(gca,'YLim',ylim)
title('Schematic')