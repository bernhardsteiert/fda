close all

load('Workspaces/site_1.mat')
timeshift = 0;
% timeshift = 125;
% timestamp = timestamp - timeshift; % Shift to main data set

% extension = '12-08-2013';
% extension = '12-25-2013';
extension = '';

sites_all = [1:70];
sites_colored = [22:24 21 2 4]; % BTC vs. IGF Native
% sites_colored = [22:24 21]; % only IGF
% sites_colored = [2:4 21]; % only BTC
% sites_colored = [16:19 21]; % only EGF

pcs = [1 2];

% f1 = figure;
% 
% xfac = 1;
% yfac = 1;
% fontsize = 20;
% 
% setFigure(f1,xfac,yfac,fontsize)
% 
% hold on
% 
% c_signal = [];
% celltype = [];
% 
% for isite = sites_all
%     [scores tmpsignal] = fPCA_late(isite,extension,timeshift);
%     c_signal = [c_signal tmpsignal];
%     celltype = [celltype ones(1,size(tmpsignal,2))*isite];
%     plot(scores(pcs(1),:),scores(pcs(2),:),'.','Color',[.7 .7 .7])
% end
% 
% ncolor = 201;
% colmap = jet(ncolor);
% dose_range = [0 100];
% dose_inds = 10.^linspace(log10(max([dose_range(1) .1])),log10(dose_range(2)),floor(ncolor/2));
% 
% firstprop = siteprop(sites_colored(1));
% for isite = sites_colored
%     scores = fPCA_late(isite,extension,timeshift);
%     sprop = siteprop(isite);
%     
%     colfac = 2*(sprop.lig_index == firstprop.lig_index)-1;
%     [tmp colind] = min(abs(sprop.lig_dose - dose_inds));
%     mycolor = colmap(ceil(ncolor/2) + colfac*colind,:);
%     
%     plot(scores(pcs(1),:),scores(pcs(2),:),'o','Color',mycolor,'MarkerFaceColor',mycolor)
% %     plotEllipsis(scores(pcs(1),:),scores(pcs(2),:),mycolor,.5);
%     
% end
% 
% axisEqual(get(gcf,'Position'))
% 
% ylabel(['PC ' num2str(pcs(2))])
% xlabel(['PC ' num2str(pcs(1))])
% 
% clim = log10(dose_range);
% clim(1) = max([clim(1) 0]);
% set(gca,'CLim',clim)
% colormap(colmap)
% colorbar('YTick',log10([1 10 100]),'YTickLabel',{sprop.lig_name,'No Stim',firstprop.lig_name}) % Vertical colorbar

% -------------------------------------------------------------------------

pcs = [3 6]; % 1-3 early; 4-6 late

f3 = figure;
hold on

xfac = 1;
yfac = 1;
fontsize = 10;

setFigure(f3,xfac,yfac,fontsize)

highdoses = [];
all_scores = [];
for isite = sites_all
    scores_early = fPCA(isite,extension,timeshift);
    scores_late = fPCA_late(isite,extension,timeshift);
    scores = [scores_early; scores_late];
%     plot(scores(pcs(1),:),scores(pcs(2),:),'.','Color',[.7 .7 .7])
    
    all_scores = [all_scores scores];
    
    sprop = siteprop(isite);
end

plot([-.3 .4],[-.3 .4],'k','LineWidth',2)
[bcoef,bint,r,rint,stats] = regress(all_scores(4,:)',[all_scores(1:3,:)' ones(size(all_scores,2),1)]);
plot(all_scores(4,:)',[all_scores(1:3,:)' ones(size(all_scores,2),1)]*bcoef,'o','MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',3)
text(0,.4,sprintf('$R^2=%1.3f$',stats(1)),'Interpreter','latex')
xlabel('Late PC1')
ylabel('$\sum_i b_i pc_i + c$','Interpreter','latex')

return

highdoses = [4 17 24 37 44 57 64];
% highdoses = [1:10]; % BTC w / wo MEKi
% highdoses = [16 17 18 19 21 12 13 14 15]; % EGF w / wo MEKi
% highdoses = [25 24 23 22 21 29 28 27 26]; % IGF w / wo MEKi

color_ind = 1;
colmap = hsv(length(highdoses));
legstr = {};
for isite = highdoses
    s = siteprop(isite);
    titstr = s.lig_name;
    titstr = sprintf('%s %g',titstr,s.lig_dose);
    if s.inh_dose > 0
%         titstr = sprintf('%s%s %i muM',titstr,s.inh_name,s.inh_dose);
        titstr = sprintf('%s%s',titstr,s.inh_name);
    end
    legstr{end+1} = titstr;
    
    scores_early = fPCA(isite,extension,timeshift);
    scores_late = fPCA_late(isite,extension,timeshift);
    scores = [scores_early; scores_late];
    plot(scores(pcs(1),:),scores(pcs(2),:),'o','Color',colmap(isite == highdoses,:),'MarkerFaceColor',colmap(isite == highdoses,:))
    plotEllipsis(scores(pcs(1),:),scores(pcs(2),:),colmap(isite == highdoses,:),.5);
end

% xlim = [-.16 .24];
% set(gca,'XLim',xlim)

% axisEqual(get(gcf,'Position'))

ylabel(['PC ' num2str(pcs(2))])
% set(gca,'YTick',-.1:.1:.2)
xlabel(['PC ' num2str(pcs(1))])
% arrow([-.13 .06],[.02 .15],'Width',.5,'Length',7)

set(gca,'CLim',[0 1])
colormap(colmap)
colorbar('YTick',linspace(1./(2*length(highdoses)),1-1./(2*length(highdoses)),length(highdoses)),'YTickLabel',legstr,'TickLength', [0 0]) % Vertical colorbar
% return
%% -------------------------------------------------------------------------

f4 = figure;

xfac = 1;
yfac = 1;
fontsize = 4;

setFigure(f4,xfac,yfac,fontsize)

rowstocols = 0.5;
nrows = 7;
ncols = 10;

% baredges = linspace(0,0.022,21); % radial_dist.m
% baredges = linspace(0,0.26,16); % edge_snr_score_pw.m dists
% baredges = linspace(0,30,16); % edge_snr_score_pw.m SNRs
% baredges = linspace(0,10,11); % edge_snr_score_pw.m nEdges
baredges = linspace(0,1.2,16); % edge_snr_score_pw_distdur.m

dists = [];
celltypeharm = [];

nEdges = [];
SNRs = [];

for isite = sites_all
    subplot(nrows,ncols,subplotpos(isite))
    
%     radial_dists = radial_dist(isite,extension,timeshift);,
    [radial_dists tmp1 tmp2 nEdge SNR] = edge_snr_score_pw_distdur(isite,extension,timeshift);
    nEdges = [nEdges nEdge];
    SNRs = [SNRs SNR];
    dists = [dists radial_dists];
    celltypeharm = [celltypeharm ones(size(radial_dists))*isite];
    
    bar(baredges,histc(dists(celltypeharm == isite),baredges));
%     bar(baredges,histc(SNRs(celltypeharm == isite),baredges));
%     bar(baredges,histc(nEdges(celltypeharm == isite),baredges));

    s = siteprop(isite);
    titstr = s.lig_name;
    titstr = sprintf('%s %g',titstr,s.lig_dose);
    if s.inh_dose > 0
%         titstr = sprintf('%s%s %i muM',titstr,s.inh_name,s.inh_dose);
        titstr = sprintf('%s%s',titstr,s.inh_name);
    end
    title(titstr)
    
    set(gca,'XLim',[-.002 1.1*max(baredges)])
end

%% -------------------------------------------------------------------------
    
return
% close all

f2 = figure;

xfac = 1;
yfac = 1;
fontsize = 4;

setFigure(f2,xfac,yfac,fontsize)

rowstocols = 0.5;
nrows = 7;
ncols = 10;

first_n = 10; % Plot first_n data-sets colored

for isite = sites_all
    subplot(nrows,ncols,subplotpos(isite))
    plot(repmat(timestamp,1,sum(celltype == isite)),c_signal(:,celltype == isite),'Color',[.7 .7 .7])
    hold on
    first_n = min(first_n,sum(celltype == isite));
    tmpind = find(celltype == isite);
    plot(repmat(timestamp,1,first_n),c_signal(:,tmpind(1:first_n)))
    plot(timestamp,nanmean(c_signal(:,celltype == isite),2),'color','k','LineWidth',2)
    
    % Comment out soon:
%     scores = fPCA(isite,extension,timeshift);
%     inds = find(celltype == isite);
%     plot(repmat(timestamp,1,length(inds(scores(3,:) < -.005))),c_signal(:,inds(scores(3,:) < -.005)))
    
    set(gca,'XLim',[50 200])
    set(gca,'YLim',[-.04 .04])
    plot([120 120],get(gca,'YLim'),'b--')
    s = siteprop(isite);
    titstr = s.lig_name;
    titstr = sprintf('%s %g',titstr,s.lig_dose);
    if s.inh_dose > 0
%         titstr = sprintf('%s%s %i muM',titstr,s.inh_name,s.inh_dose);
        titstr = sprintf('%s%s',titstr,s.inh_name);
    end
    title(titstr)
end

%% -------------------------------------------------------------------------
    
% return
% close all

f5 = figure;

xfac = 1;
yfac = 1;
fontsize = 4;

setFigure(f5,xfac,yfac,fontsize)

rowstocols = 0.5;
nrows = 7;
ncols = 10;

ncolor = 201;
color = repmat(linspace(0,1,ncolor),3,1);

color = color(:,end:-1:1); % Gray scale - darkness depending on score

[radial_dist_sorted ind_sort_radial] = sort(dists);
radial_space = linspace(min(dists),dists(ind_sort_radial(ceil(length(dists)*.99))),ncolor);

for isite = sites_all
    subplot(nrows,ncols,subplotpos(isite))
    hold on
    
    c_signal_single = c_signal(:,ind_sort_radial);
    c_signal_single = c_signal_single(:,(celltypeharm(ind_sort_radial) == isite));
    radial_dist_single = radial_dist_sorted(celltypeharm(ind_sort_radial) == isite);
    
    for i = 1:size(c_signal_single,2)
        [tmp color_ind] = min(abs(radial_dist_single(i) - radial_space));
        plot(timestamp,c_signal_single(:,i),'Color',color(:,color_ind))
    end
    
    plot(get(gca,'XLim'),[0 0],'--k')
    
    set(gca,'XLim',[50 510])
    
    plot([120 120],[-0.04 0.04],'b--')
    
    if subplotpos(isite) == (nrows-1)*ncols+1
        xlabel('time [min]')
    end
    if subplotpos(isite) == 2
        ylabel('log_{10} FOXO3a Cyt/Nuc ratio');
    end
    set(gca,'YLim',[-0.04 0.04])
    s = siteprop(isite);
    titstr = s.lig_name;
    titstr = sprintf('%s %g',titstr,s.lig_dose);
    if s.inh_dose > 0
%         titstr = sprintf('%s%s %i muM',titstr,s.inh_name,s.inh_dose);
        titstr = sprintf('%s%s',titstr,s.inh_name);
    end
    title(titstr)
end