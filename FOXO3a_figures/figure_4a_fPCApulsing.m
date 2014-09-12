addpath('./Functions/')
load('./Workspaces/scores_04-15_new')

extension = '04-15-2014';
pcs = [2 3];

highdoses_all = [1 24 25 48 49 66; 3 22 27 46 51 65; 5 20 29 44 53 64; 7 18 31 42 55 63; 9 16 33 40 57 62; 11 14 35 38 59 61];

kswidth_all = ones(1,6)*35;
puls_thres = [.3 .5 .45 .2 .2 .2];

for icell = 1:size(highdoses_all,1)

    highdoses = highdoses_all(icell,:);
    kswidth = kswidth_all(icell);

    highinds = ismember(celltype,highdoses);
    [dists_boxcox,lambda,c] = boxcox(dists(highinds));
    thres_trafo = boxcox_apply(puls_thres(icell),lambda,c);
    mymin = min(dists_boxcox);
    myrange = range(dists_boxcox);
    dists_boxcox = (dists_boxcox-mymin) ./ myrange;
    thres_trafo = (thres_trafo-mymin) ./ myrange;

    figure

    plot(scores_all(pcs(1),:),scores_all(pcs(2),:),'.','Color',[.7 .7 .7])
    hold on

    color_ind = 1;
    colmap = hsv(length(highdoses));
    legstr = {};
    for isite = highdoses
        s = siteprop(isite,extension);
        titstr = s.lig_name;
    %     titstr = sprintf('%s + %s %g + %s %g - %s',titstr,s.drug1_name,s.drug1_dose,s.drug2_name,s.drug2_dose,s.celltype);
        legstr{end+1} = titstr;

        plot(scores_all(pcs(1),celltype == isite),scores_all(pcs(2),celltype == isite),'o','Color',colmap(isite == highdoses,:),'MarkerFaceColor',colmap(isite == highdoses,:))
        plotEllipsis(scores_all(pcs(1),celltype == isite),scores_all(pcs(2),celltype == isite),colmap(isite == highdoses,:),.5);
    end

    xlim = [-.4 .3];
    set(gca,'XLim',xlim)
    ylim = [-.08 .12];
    set(gca,'YLim',ylim)

    % axisEqual(get(gcf,'Position'))

    title(s.celltype)
    ylabel(['PC ' num2str(pcs(2))])
    % set(gca,'YTick',-.1:.1:.2)
    xlabel(['PC ' num2str(pcs(1))])

    set(gca,'CLim',[0 1])
    colormap(colmap)
    colorbar('YTick',linspace(1./(2*length(highdoses)),1-1./(2*length(highdoses)),length(highdoses)),'YTickLabel',legstr,'TickLength', [0 0]) % Vertical colorbar


    figure
    hold on

    legh = [];
    legstr = {};
    for isite = highdoses
        s = siteprop(isite,extension);

        [f,xi] = ksdensity(dists_boxcox(celltype(highinds) == isite),'width',range(dists)./kswidth);
        legh = [legh plot(xi,f,'Color',colmap(isite == highdoses,:))];

        legstr{end+1} = s.lig_name;

    end
    % set(gca,'XLim',[min(dists_boxcox) max(dists_boxcox)]+[-.2 .2]*range(dists_boxcox))
    set(gca,'XLim',[-.2 1.2])
    plot([thres_trafo thres_trafo],get(gca,'YLim'),'k--')
    title(s.celltype)
    legend(legh,legstr)

    xlabel('pulsatory score')
    ylabel('probability density')
    
end