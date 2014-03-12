extension = '01-22-2014';

% Initialize timestamp
remotepath = mypath();
    
warning('off','MATLAB:dispatcher:pathWarning')

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
addpath(grabdataPath)

if exist(remotepath,'dir')
    timestamp = grabdata_new(isite,extension);
end

close all

% load('Workspaces/site_1_01-22-2014.mat')
timeshift = -70;
timestamp = timestamp - timeshift; % Shift to main data set

% sites_all = [1:57 59 60];
sites_all = [1:12 49:57 59 60];
% sites_colored = [49:54]; % No Boretezo
sites_colored = [1:12];


pcs = [2 3];

% f1 = figure;
% 
xfac = 1;
yfac = 1;
fontsize = 14;
% 
% setFigure(f1,xfac,yfac,fontsize)
% 
% hold on
% 
% c_signal = [];
% celltype = [];
% 
% for isite = sites_all
%     [scores tmpsignal] = fPCA(isite,extension,timeshift);
%     c_signal = [c_signal tmpsignal];
%     celltype = [celltype ones(1,size(tmpsignal,2))*isite];
%     plot(scores(pcs(1),:),scores(pcs(2),:),'.','Color',[.7 .7 .7])
% end
% 
% % return
% 
% ncolor = 201;
% colmap = jet(ncolor);
% dose_range = [0 100];
% dose_inds = 10.^linspace(log10(max([dose_range(1) .1])),log10(dose_range(2)),floor(ncolor/2));
% 
% firstprop = siteprop(sites_colored(1),extension);
% for isite = sites_colored
%     scores = fPCA(isite,extension,timeshift);
%     sprop = siteprop(isite,extension);
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
% set(gca,'YLim',[-.05 .03])
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

f3 = figure;
hold on

xfac = 1;
yfac = 1;

setFigure(f3,xfac,yfac,fontsize)

highdoses = [];
c_signal = [];
celltype = [];

for isite = sites_all
    [scores tmpsignal] = fPCA(isite,extension,timeshift);
    c_signal = [c_signal tmpsignal];
    celltype = [celltype ones(1,size(tmpsignal,2))*isite];
    plot(scores(pcs(1),:),scores(pcs(2),:),'.','Color',[.7 .7 .7])
end


% highdoses = [4 17 24 21]; % BTC EGF IGF NS
% highdoses = [49:54]; % No Boretezo
% highdoses = [5 1 3 4 6 2];
highdoses = 7:12;

color_ind = 1;
colmap = hsv(length(highdoses));
legstr = {};
for isite = highdoses
    s = siteprop(isite,extension);
    titstr = [s.celltype ': ' s.lig_name];
    titstr = sprintf('%s %g',titstr,s.lig_dose);
    titstr = [titstr s.drug_name ' ' s.drug_timing];
    legstr{end+1} = titstr;
    
    scores = fPCA(isite,extension,timeshift);
    plot(scores(2,:),scores(3,:),'o','Color',colmap(isite == highdoses,:),'MarkerFaceColor',colmap(isite == highdoses,:))
    plotEllipsis(scores(2,:),scores(3,:),colmap(isite == highdoses,:),.5);
end

xlim = [-.5 .3];
set(gca,'XLim',xlim)
set(gca,'YLim',[-.2 .3])

% axisEqual(get(gcf,'Position'))

ylabel(['PC ' num2str(pcs(2))])
% set(gca,'YTick',-.1:.1:.2)
xlabel(['PC ' num2str(pcs(1))])
% arrow([-.13 .06],[.02 .15],'Width',.5,'Length',7)

set(gca,'CLim',[0 1])
colormap(colmap)
colorbar('YTick',linspace(1./(2*length(highdoses)),1-1./(2*length(highdoses)),length(highdoses)),'YTickLabel',legstr,'TickLength', [0 0]) % Vertical colorbar
% return
% -------------------------------------------------------------------------

f4 = figure;

xfac = 1;
yfac = 1;
fontsize = 4;

setFigure(f4,xfac,yfac,fontsize)

rowstocols = 0.5;
nrows = 5;
ncols = 12;

% baredges = linspace(0,0.022,21); % radial_dist.m
% baredges = linspace(0,0.12,16); % edge_snr_score_pw.m
baredges = linspace(0,2,16);

dists = [];
celltypeharm = [];

for isite = sites_all
    subplot(nrows,ncols,subplotpos(isite,ncols))
    
%     radial_dists = radial_dist(isite,extension,timeshift);,
%     radial_dists = edge_snr_score_pw(isite,extension,timeshift);
    radial_dists = edge_snr_score_pw_distdur(isite,extension,timeshift,1/60);
    dists = [dists radial_dists];
    celltypeharm = [celltypeharm ones(size(radial_dists))*isite];
    
    bar(baredges,histc(dists(celltypeharm == isite),baredges));

    s = siteprop(isite,extension);
    titstr = [s.celltype '\newline' s.lig_name];
    if s.lig_dose > 0
        titstr = sprintf('%s%g',titstr,s.lig_dose);
    end
    if ~isempty(s.drug_timing)
        titstr = [titstr ' ' s.drug_name(1:2) ' ' s.drug_timing(1:7)];
    end
    title(titstr)
    
    set(gca,'XLim',[-.002 1.1*max(baredges)])
end

% -------------------------------------------------------------------------
    
% return
% close all

f2 = figure;

xfac = 1;
yfac = 1;
fontsize = 4;

setFigure(f2,xfac,yfac,fontsize)

rowstocols = 0.5;
nrows = 5;
ncols = 12;

first_n = 10; % Plot first_n data-sets colored

for isite = sites_all
    subplot(nrows,ncols,subplotpos(isite,ncols))
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
%     set(gca,'YLim',[-.01 .01])
    plot([120 120],get(gca,'YLim'),'b--')
    s = siteprop(isite,extension);
    titstr = [s.celltype '\newline' s.lig_name];
    if s.lig_dose > 0
        titstr = sprintf('%s%g',titstr,s.lig_dose);
    end
    if ~isempty(s.drug_timing)
        titstr = [titstr ' ' s.drug_name(1:2) ' ' s.drug_timing(1:7)];
    end
    title(titstr)
end

% -------------------------------------------------------------------------
    
% return
% close all

f5 = figure;

xfac = 1;
yfac = 1;
fontsize = 4;

setFigure(f5,xfac,yfac,fontsize)

rowstocols = 0.5;
nrows = 5;
ncols = 12;

ncolor = 201;
color = repmat(linspace(0,1,ncolor),3,1);

color = color(:,end:-1:1); % Gray scale - darkness depending on score

[radial_dist_sorted ind_sort_radial] = sort(dists);
radial_space = linspace(min(dists),dists(ind_sort_radial(ceil(length(dists)*.99))),ncolor);

for isite = sites_all
    subplot(nrows,ncols,subplotpos(isite,ncols))
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
    
    plot([120 120],[-0.06 0.04],'b--')
    
    if subplotpos(isite) == (nrows-1)*ncols+1
        xlabel('time [min]')
    end
    if subplotpos(isite) == 2
        ylabel('log_{10} FOXO3a Cyt/Nuc ratio');
    end
%     set(gca,'YLim',[-0.06 0.04])
    s = siteprop(isite,extension);
    titstr = [s.celltype '\newline' s.lig_name];
    if s.lig_dose > 0
        titstr = sprintf('%s%g',titstr,s.lig_dose);
    end
    if ~isempty(s.drug_timing)
        titstr = [titstr ' ' s.drug_name(1:2) ' ' s.drug_timing(1:7)];
    end
    title(titstr)
end