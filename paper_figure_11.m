close all

icount = 1;

for idose = [2.5 5 10 20 50 100]

    f2 = subplot(2,3,icount);
    icount = icount + 1;
    hold on

    xfac = 1;
    yfac = .6;

%     setFigure(f2,xfac,yfac,fontsize)
    
    title(num2str(idose,3))

    highdoses = [];
    for isite = sites_all
        scores = fPCA(isite);
        plot(scores(2,:),scores(3,:),'.','Color',[.7 .7 .7])

        sprop = siteprop(isite);
        if sprop.lig_dose == idose
            highdoses = [highdoses isite];
        end
    end

    resort = [2 3 4 1 6 5];
    highdoses = highdoses(resort);

    color_ind = 1;
    colmap = hsv(length(highdoses));
    legstr = {};
    for isite = highdoses
        sprop = siteprop(isite);
        legstr{end+1} = sprop.lig_name(1:3);

        scores = fPCA(isite);
        plot(scores(2,:),scores(3,:),'.','Color',colmap(isite == highdoses,:))
        plotEllipsis(scores(2,:),scores(3,:),colmap(isite == highdoses,:),.5);
    end

    xlim = [-.16 .24];
    set(gca,'XLim',xlim)
    set(gca,'YLim',[-.09 .13])

%     axisEqual(get(f2,'Position'))

    ylabel(['PC ' num2str(pcs(2))])
    % arrow([-.13 -.1],[.2 -.1],'Width',.5,'Length',7)
    set(gca,'YTick',-.1:.1:.2)
    xlabel(['PC ' num2str(pcs(1))])
    % arrow([-.13 .06],[.02 .15],'Width',.5,'Length',7)

    set(gca,'CLim',[0 1])
    colormap(colmap)
%     colorbar('YTick',linspace(1./(2*length(highdoses)),1-1./(2*length(highdoses)),length(highdoses)),'YTickLabel',legstr,'TickLength', [0 0]) % Vertical colorbar

end