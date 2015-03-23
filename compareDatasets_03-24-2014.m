close all

load('Workspaces/site_1_03-24-2014.mat')
timeshift = 52;
timestamp = timestamp - timeshift; % Shift to main data set

extension = '03-24-2014';

sites_all = [1:60];
% sites_all = [1:10 12:34 36:58 62:66];
sites_colored = [4 17 24 37 44 57]; % BTC vs. IGF Native
% sites_colored = [22:24 21]; % only IGF
% sites_colored = [2:4 21]; % only BTC
% sites_colored = [16:19 21]; % only EGF

pcs = [2 3];

f1 = figure;

xfac = 1;
yfac = 1;
fontsize = 20;

setFigure(f1,xfac,yfac,fontsize)

hold on

c_signal = [];
celltype = [];

for isite = sites_all
    [scores tmpsignal] = fPCA_corrected(isite,extension,timeshift);
    c_signal = [c_signal tmpsignal];
    celltype = [celltype ones(1,size(tmpsignal,2))*isite];
    plot(scores(pcs(1),:),scores(pcs(2),:),'.','Color',[.7 .7 .7])
end

c_signal(isinf(c_signal)) = nan;

ncolor = 201;
colmap = jet(ncolor);
dose_range = [0 100];
dose_inds = 10.^linspace(log10(max([dose_range(1) .1])),log10(dose_range(2)),floor(ncolor/2));

firstprop = siteprop(sites_colored(1));
for isite = sites_colored
    scores = fPCA_corrected(isite,extension,timeshift);
    sprop = siteprop(isite);
    
    colfac = 2*(sprop.lig_index == firstprop.lig_index)-1;
    [tmp colind] = min(abs(sprop.lig_dose - dose_inds));
    mycolor = colmap(ceil(ncolor/2) + colfac*colind,:);
    
    plot(scores(pcs(1),:),scores(pcs(2),:),'o','Color',mycolor,'MarkerFaceColor',mycolor)
%     plotEllipsis(scores(pcs(1),:),scores(pcs(2),:),mycolor,.5);
    
end

axisEqual(get(gcf,'Position'))

ylabel(['PC ' num2str(pcs(2))])
xlabel(['PC ' num2str(pcs(1))])

clim = log10(dose_range);
clim(1) = max([clim(1) 0]);
set(gca,'CLim',clim)
colormap(colmap)
colorbar('YTick',log10([1 10 100]),'YTickLabel',{sprop.lig_name,'No Stim',firstprop.lig_name}) % Vertical colorbar

% -------------------------------------------------------------------------

f3 = figure;
hold on

xfac = 1;
yfac = 1;

setFigure(f3,xfac,yfac,fontsize)

highdoses = [];
scores_all = [];
for isite = sites_all
    scores = fPCA_corrected(isite,extension,timeshift);
    scores_all = [scores_all scores];
end

plot(scores_all(2,:),scores_all(3,:),'.','Color',[.7 .7 .7])

outliers = logical(zeros(size(celltype)));
outlier_thres = .05;
nNeighbours = 8;
for iscore = 1:length(scores_all)
    sorted_dists = sort(sqrt(sum((repmat(scores_all(2:3,iscore),1,size(scores_all,2))-scores_all(2:3,:)).^2,1)));
    outliers(iscore) = sorted_dists(nNeighbours+1) > outlier_thres;
end

plot(scores_all(2,outliers),scores_all(3,outliers),'ro')
%
f3 = figure;
hold on

xfac = 1;
yfac = 1;

setFigure(f3,xfac,yfac,fontsize)

highdoses = [11 14 35 38 59 60]; % 184A1ERKspec

plot(scores_all(2,~outliers),scores_all(3,~outliers),'.','Color',[.7 .7 .7])

color_ind = 1;
colmap = hsv(length(highdoses));
legstr = {};
for isite = highdoses
    s = siteprop(isite,extension);
    titstr = s.lig_name;
%     titstr = sprintf('%s + %s %g + %s %g - %s',titstr,s.drug1_name,s.drug1_dose,s.drug2_name,s.drug2_dose,s.celltype);
    legstr{end+1} = titstr;
    
    scores = fPCA_corrected(isite,extension,timeshift);
    plot(scores(2,~outliers(celltype == isite)),scores(3,~outliers(celltype == isite)),'o','Color',colmap(isite == highdoses,:),'MarkerFaceColor',colmap(isite == highdoses,:))
    plotEllipsis(scores(2,~outliers(celltype == isite)),scores(3,~outliers(celltype == isite)),colmap(isite == highdoses,:),.5);
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

w = 10;
c_signal_slidwin = nan(size(c_signal,1)-2*w+1,size(c_signal,2));
for it = (w+1):size(c_signal,1)-w
    c_signal_slidwin(it-w,:) = c_signal(it,:)-nanmean(c_signal((it-w):(it+w),:),1);
end

c_signal = c_signal_slidwin;
timestamp = timestamp(1:(end-2*w+1));


%%

f4 = figure;

xfac = 1;
yfac = 1;
fontsize = 4;

setFigure(f4,xfac,yfac,fontsize)

rowstocols = 0.5;
nrows = 6;
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
amps = [];
peakdurs = [];

for isite = sites_all
    subplot(nrows,ncols,subplotpos(isite,ncols))
    
%     radial_dists = radial_dist(isite,extension,timeshift);,
    [radial_dists tmp1 tmp2 nEdge SNR amp tmp3 peakdur] = edge_snr_score_pw_distdur(isite,extension,timeshift,1/120,'harm_basis_04-15-2014_late',1,1);
    nEdges = [nEdges nEdge];
    SNRs = [SNRs SNR];
    dists = [dists radial_dists];
    amps = [amps amp];
    peakdurs = [peakdurs peakdur];
    celltypeharm = [celltypeharm ones(size(radial_dists))*isite];
    
    bar(baredges,histc(dists(celltypeharm == isite),baredges));
%     bar(baredges,histc(SNRs(celltypeharm == isite),baredges));
%     bar(baredges,histc(nEdges(celltypeharm == isite),baredges));

    s = siteprop(isite,extension);
    titstr = s.lig_name;
%     titstr = sprintf('%s + %s %g + %s %g - %s',titstr,s.drug1_name,s.drug1_dose,s.drug2_name,s.drug2_dose,s.celltype);
    title(titstr)
    
    set(gca,'XLim',[-.002 1.1*max(baredges)])
end

%% -------------------------------------------------------------------------
    
% return 
% close all

f2 = figure;

xfac = 1;
yfac = 1;
fontsize = 4;

setFigure(f2,xfac,yfac,fontsize)

rowstocols = 0.5;
nrows = 6;
ncols = 10;

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
%     scores = fPCA_corrected(isite,extension,timeshift);
%     inds = find(celltype == isite);
%     plot(repmat(timestamp,1,length(inds(scores(3,:) < -.005))),c_signal(:,inds(scores(3,:) < -.005)))
    
    set(gca,'XLim',[50 500])
    set(gca,'YLim',[-.2 .2])
    plot([120 120],get(gca,'YLim'),'b--')
    s = siteprop(isite,extension);
    titstr = s.lig_name;
%     titstr = sprintf('%s + %s %g + %s %g - %s',titstr,s.drug1_name,s.drug1_dose,s.drug2_name,s.drug2_dose,s.celltype);
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
nrows = 6;
ncols = 10;

ncolor = 201;
color = repmat(linspace(0,1,ncolor),3,1);

color = color(:,end:-1:1); % Gray scale - darkness depending on score

[radial_dist_sorted ind_sort_radial] = sort(dists);
radial_space = linspace(min(dists),dists(ind_sort_radial(ceil(sum(~isnan(dists))*.99))),ncolor);

for isite = sites_all
    subplot(nrows,ncols,subplotpos(isite,ncols))
    hold on
    
    c_signal_single = c_signal(:,ind_sort_radial);
    c_signal_single = c_signal_single(:,(celltypeharm(ind_sort_radial) == isite));
    radial_dist_single = radial_dist_sorted(celltypeharm(ind_sort_radial) == isite);
    
%     noNaN = sum(isnan(c_signal_single),1)==0;
%     c_signal_single = c_signal_single(:,noNaN);
%     radial_dist_single = radial_dist_single(noNaN);
    
    for i = 1:size(c_signal_single,2)
        [tmp color_ind] = min(abs(radial_dist_single(i) - radial_space));
        plot(timestamp,c_signal_single(:,i),'Color',color(:,color_ind))
    end
    
    plot(get(gca,'XLim'),[0 0],'--k')
    
    set(gca,'XLim',[50 800])
    
    plot([120 120],[-0.2 0.2],'b--')
    
    if subplotpos(isite) == (nrows-1)*ncols+1
        xlabel('time [min]')
    end
    if subplotpos(isite) == 2
        ylabel('log_{10} FOXO3a Cyt/Nuc ratio');
    end
    set(gca,'YLim',[-0.2 0.2])
    s = siteprop(isite,extension);
    titstr = s.lig_name;
%     titstr = sprintf('%s + %s %g + %s %g - %s',titstr,s.drug1_name,s.drug1_dose,s.drug2_name,s.drug2_dose,s.celltype);
    title(titstr)
end

%% Pulsing(Ligand dose)
% Run cell below --> ratio_puls_mat
close all
imagesc(ratio_puls_mat')
colorbar
possible_doses = [100 100/2.5 100/2.5^2 100/2.5^3 100/2.5^4 100/2.5^5 100/2.5^6 100/2.5^7 100/2.5^8 0];
possible_ligands = {'IGF','HRG','HGF','EGF','BTC','EPR'};
set(gca,'XTick',1:length(possible_doses),'XTickLabel',possible_doses)
set(gca,'YTick',1:length(possible_ligands),'YTickLabel',possible_ligands)

possible_doses_main = [100 50 20 10 5 2.5 0];
log_doses = log10(possible_doses);
log_doses(end) = log10(100/2.5^9);
log_doses_main = log10(possible_doses_main);
log_doses_main(end) = log_doses(end);

ratio_puls_mat_interp = nan(size(ratio_puls_mat,2),length(log_doses_main));
for ilig = 1:size(ratio_puls_mat,2)
    ratio_puls_mat_interp(ilig,:) = interp1(log_doses,ratio_puls_mat(:,ilig)',log_doses_main);
end

figure
imagesc(ratio_puls_mat_interp)
colorbar
set(gca,'XTick',1:length(possible_doses_main),'XTickLabel',possible_doses_main)
set(gca,'YTick',1:length(possible_ligands),'YTickLabel',possible_ligands)

%% Sort by ratio of pulsing cells
f6 = figure;

xfac = 1;
yfac = 1;
fontsize = 16;

setFigure(f6,xfac,yfac,fontsize)

puls_thres = .8;

labels = {};
ratio_puls = nan(size(sites_all));
ratio_puls_mat = nan(ncols,nrows);

mean_amp = nan(size(sites_all));
mean_amp_mat = nan(ncols,nrows);

mean_peakdur = nan(size(sites_all));
mean_peakdur_mat = nan(ncols,nrows);

for isite = sites_all
    ratio_puls(isite) = sum(dists(celltype == isite) > puls_thres) / sum(~isnan(dists(celltype == isite)));
    mean_amp(isite) = mean(amps(dists > puls_thres & celltype == isite));
    mean_peakdur(isite) = mean(peakdurs(dists > puls_thres & celltype == isite));
%     mean_amp(isite) = nanmean(amps(celltype == isite));
%     mean_peakdur(isite) = nanmean(peakdurs(celltype == isite));
    ratio_puls_mat(subplotpos(isite,ncols)) = ratio_puls(isite);
    mean_amp_mat(subplotpos(isite,ncols)) = mean_amp(isite);
    mean_peakdur_mat(subplotpos(isite,ncols)) = mean_peakdur(isite);
    
    s = siteprop(isite,extension);
    titstr = s.lig_name;
%     titstr = sprintf('%s %g + %s %g - %s',titstr,s.lig_dose,s.drug_name,s.drug_dose,s.celltype);
    labels{end+1} = titstr;
end

[ratio_puls_sorted ind_puls_sorted] = sort(ratio_puls,'descend');

[X,Y] = meshgrid(1:size(ratio_puls_mat,2), 1:1:size(ratio_puls_mat,1)); 

valid = ~isnan(ratio_puls_mat); 
M = griddata(X(valid),Y(valid),ratio_puls_mat(valid),X,Y);
zlab = 'Ratio of pulsing cells';

% valid = ~isnan(mean_amp_mat); 
% M = griddata(X(valid),Y(valid),mean_amp_mat(valid),X,Y); 
% zlab = 'Amplitude';

% valid = ~isnan(mean_peakdur_mat); 
% M = griddata(X(valid),Y(valid),mean_peakdur_mat(valid),X,Y); 
% zlab = 'Peak duration';

% surf(M(1:6,:))
imagesc(M(1:6,:))
set(gca,'XDir','Reverse')
title('MCF10A')
xlabel('EGF dose')
ylabel('MEKi dose')
zlabel(zlab)
set(gca,'XTick',1:6,'XTickLabel',[100 20 4 .8 .16 0])
set(gca,'YTick',1:6,'YTickLabel',[.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0])

f8 = figure;

xfac = 1;
yfac = 1;
fontsize = 16;

setFigure(f8,xfac,yfac,fontsize)

% surf(M(12:-1:7,:))
imagesc(M(12:-1:7,:))
set(gca,'XDir','Reverse')
title('184A1')
xlabel('EGF dose')
ylabel('MEKi dose')
zlabel(zlab)
set(gca,'XTick',1:6,'XTickLabel',[100 20 4 .8 .16 0])
set(gca,'YTick',1:6,'YTickLabel',[.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0])

% bar(ratio_puls_sorted)
% set(gca,'XLim',[0 73])
% 
% h = my_xticklabels(sites_all,labels(ind_puls_sorted), ...
%       'Rotation',-90, ...
%       'VerticalAlignment','middle', ...
%       'HorizontalAlignment','left');
%   
% hold on
% plot(get(gca,'XLim'),[.2 .2],'k--')
% ylabel('fraction of pulsing cells')

%% Plot single traces sorted by pulsing score for given condition
% mycond = 67; % 67 = EGF 0; MEKi 0; MCF10A
% mycond = 6; % 6 = EGF 100; MEKi 0; MCF10A
mycond = 12; % 6 = EGF 100; MEKi 0.1; 184A1

nrows = 10;
ncols = 10;

myind = find(celltype == mycond);
myind = myind(1:min([100 length(myind)]));

[tmp indsorted] = sort(dists(myind),'descend');

f8 = figure;

xfac = 1;
yfac = 1;
fontsize = 4;

setFigure(f8,xfac,yfac,fontsize)

for isite = 1:length(myind)
    subplot(nrows,ncols,isite)
    plot(timestamp,c_signal(:,myind(indsorted(isite))))
    title(sprintf('score = %g; ind = %i',dists(myind(indsorted(isite))),myind(indsorted(isite))-myind(1)+1))
    set(gca,'XLim',[0 1000])
end


%% Pulsing vs. mean response plot
load('harm_basis_corrected');
basis_eval = eval_basis(linspace(51,120,101),harm_basis);
mean_early = mean(basis_eval(:,1));

load('harm_basis_04-18-2014_late');
basis_range = getbasisrange(harm_basis);
basis_eval = eval_basis(linspace(basis_range(1),basis_range(2),101),harm_basis);
mean_late = mean(basis_eval(:,1));

mean_response = [];
for isite = sites_all
    
    scores_early = fPCA_corrected(isite,extension,timeshift,'harm_basis_corrected',1);
    scores_late = fPCA_corrected(isite,extension,timeshift,'harm_basis_04-18-2014_late',1);
    mean_response = [mean_response scores_late(1,:)*mean_late-scores_early(1,:)*mean_early];
    
end

%%
f6 = figure;

xfac = 1;
yfac = 1;
fontsize = 16;

setFigure(f6,xfac,yfac,fontsize)

labels = {};
mean_response_mat = nan(ncols,nrows);

for isite = sites_all
    mean_response_mat(isite) = mean(mean_response(celltype==isite));
end


[X,Y] = meshgrid(1:size(mean_response_mat,2), 1:1:size(mean_response_mat,1)); 

valid = ~isnan(mean_response_mat); 
M = griddata(X(valid),Y(valid),mean_response_mat(valid),X,Y); 

surf(M(1:6,:))
title('MCF10A')
xlabel('EGF dose')
ylabel('MEKi dose')
set(gca,'XTick',1:6,'XTickLabel',[100 20 4 .8 .16 0])
set(gca,'YTick',1:6,'YTickLabel',[.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0])
view(-45,60)

f8 = figure;

xfac = 1;
yfac = 1;
fontsize = 16;

setFigure(f8,xfac,yfac,fontsize)

surf(M(12:-1:7,:))
title('184A1')
xlabel('EGF dose')
ylabel('MEKi dose')
set(gca,'XTick',1:6,'XTickLabel',[100 20 4 .8 .16 0])
set(gca,'YTick',1:6,'YTickLabel',[.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0])
view(-45,60)

%%
f7 = figure;

xfac = 1;
yfac = 1;
fontsize = 16;

setFigure(f7,xfac,yfac,fontsize)
hold on

for isite = sites_all
    
    s = siteprop(isite,extension);
    if strmatch(s.celltype,'MCF10A','exact')
        marker = 's';
        col1 = 1/3;
    else
        marker = 'o';
        col1 = 1;
    end
    
        
    switch s.lig_dose
        case 100
            markersize = 16;
        case 20
            markersize = 14;
        case 4
            markersize = 12;
        case .8
            markersize = 10;
        case .16
            markersize = 8;
        case 0
            markersize = 6;
    end
    
    switch s.drug_dose
        case 0
            col3 = 1/6;
        case .1/4^4
            col3 = 2/6;
        case .1/4^3
            col3 = 3/6;
        case .1/4^2
            col3 = 4/6;
        case .1/4
            col3 = 5/6;
        case .1
            col3 = 1;
    end
    
    color = hsv2rgb([col1 1 col3]);
    plot(mean(mean_response(celltype == isite)),sum(dists(celltype == isite) > puls_thres) / sum(celltype == isite),marker,'MarkerEdgeColor','k','MarkerSize',markersize,'MarkerFaceColor',color)
    
    
end

set(gca,'YLim',[-.05 .5])
xlabel('mean response (late PC1 - early PC1)')
ylabel('ratio pulsing cells')

%% Mean responses for (fixed EGF; varying MEKi) or other way round
figure

smooth_level = 500; % the larger the more smoothing
mean_signal = nan(length(timestamp),sites_all(end));
std_signal = nan(length(timestamp),sites_all(end));
n_signal = nan(length(timestamp),sites_all(end));
for isite = sites_all
    mean_signal(:,isite) = csaps(timestamp,nanmean(c_signal(:,celltype == isite),2),1/smooth_level,timestamp);
    std_signal(:,isite) = csaps(timestamp,nanstd(c_signal(:,celltype == isite),[],2),1/smooth_level,timestamp);
    n_signal(isite) = length(c_signal(:,celltype == isite));
end

fixedEGF = 13:24;
fixedMEKi = [4 24-4+1 24+4 48-4+1 48+4 72-4+1 9 24-9+1 24+9 48-9+1 48+9 72-9+1];
fixedEGF = fixedEGF([6:-1:1 7:12]);
fixedMEKi = fixedMEKi([6:-1:1 7:12]);

ligand_dose = [100 20 4 .8 .16 0];
drug_dose = [.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0];

colmap = jet(length(drug_dose));
linewidth = 1;
legh = [];
legstr = {};
for isite = fixedEGF
    sprop = siteprop(isite,extension);
    if strmatch(sprop.celltype,'MCF10A','exact')
        subplot(2,2,1)
        hold on
        legh = [legh plot(timestamp,mean_signal(:,isite),'Color',colmap(sprop.drug_dose == drug_dose,:))];
        legstr{end+1} = sprintf('MEKi %g',drug_dose(sprop.drug_dose == drug_dose));
    else
        subplot(2,2,2)
        hold on
        plot(timestamp,mean_signal(:,isite),'Color',colmap(sprop.drug_dose == drug_dose,:))
    end
    
    tmpx = [timestamp; flipud(timestamp)];
    tmpy = [mean_signal(:,isite) + 2*std_signal(:,isite)/sqrt(n_signal(isite)); flipud(mean_signal(:,isite) - 2*std_signal(:,isite)/sqrt(n_signal(isite)))];
    
    ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp, 'FaceColor', colmap(sprop.drug_dose == drug_dose,:)*0.1+0.9, 'EdgeColor', colmap(sprop.drug_dose == drug_dose,:)*0.1+0.9,'LineWidth',linewidth);
    ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', colmap(sprop.drug_dose == drug_dose,:)*0.3+0.7,'LineWidth',linewidth);
    
end
legend(legh,legstr)

legh = [];
legstr = {};
colmap = jet(length(ligand_dose));
for isite = fixedMEKi
    sprop = siteprop(isite,extension);
    if strmatch(sprop.celltype,'MCF10A','exact')
        subplot(2,2,3)
        hold on
        legh = [legh plot(timestamp,mean_signal(:,isite),'Color',colmap(sprop.lig_dose == ligand_dose,:))];
        legstr{end+1} = sprintf('EGF %g',ligand_dose(sprop.lig_dose == ligand_dose));
    else
        subplot(2,2,4)
        hold on
        plot(timestamp,mean_signal(:,isite),'Color',colmap(sprop.lig_dose == ligand_dose,:))
    end
    
    tmpx = [timestamp; flipud(timestamp)];
    tmpy = [mean_signal(:,isite) + 2*std_signal(:,isite)/sqrt(n_signal(isite)); flipud(mean_signal(:,isite) - 2*std_signal(:,isite)/sqrt(n_signal(isite)))];
    
    ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp, 'FaceColor', colmap(sprop.lig_dose == ligand_dose,:)*0.1+0.9, 'EdgeColor', colmap(sprop.lig_dose == ligand_dose,:)*0.1+0.9,'LineWidth',linewidth);
    ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', colmap(sprop.lig_dose == ligand_dose,:)*0.3+0.7,'LineWidth',linewidth);
    
end
legend(legh,legstr)

subplot(2,2,1)
title('MCF10A - EGF 20; varying MEKi')
set(gca,'YLim',[-.012 .022])
set(gca,'XLim',[0 1000])
ylabel('log_{10} FOXO3a [Cyt/Nuc]');
subplot(2,2,2)
title('184A1 - EGF 20; varying MEKi')
set(gca,'YLim',[-.012 .022])
set(gca,'XLim',[0 1000])
subplot(2,2,3)
title('MCF10A - MEKi 1/640; varying EGF')
set(gca,'YLim',[-.012 .022])
set(gca,'XLim',[0 1000])
xlabel('time [min]')
subplot(2,2,4)
title('184A1 - MEKi 1/640; varying EGF')
set(gca,'YLim',[-.012 .022])
set(gca,'XLim',[0 1000])

%% c_signal smoothing
nbasis = round(length(timestamp)/3);
basis = create_bspline_basis([timestamp(1) timestamp(end)], nbasis);
smoothed_data_woNharm = smooth_basis(timestamp,c_signal,basis);
c_signal = eval_fd(smoothed_data_woNharm,timestamp);

%% Outlier removal (not really usefull :()
for i = 1:size(c_signal,2)

    mysignal = c_signal(:,i);
    testdiff = diff(mysignal);
    testdiff(isinf(testdiff)) = nan;


    win = 15;

    allranges = nan(1,length(testdiff)-win);
    for islide = 1:length(testdiff)-win
        allranges(islide) = range(testdiff(islide:islide+win-1));
    end

    range_sorted = sort(allranges);
    range_thres = range_sorted(round(length(range_sorted)*.5))-min(range_sorted);

    outliers = find(allranges > min(range_sorted)+4*range_thres);
    cleaned_inds = ones(1,length(mysignal));

    if ~isempty(outliers)
        outlier_bd = find(diff(outliers) > 1);

        outlier_lb = [outliers(1) outliers(outlier_bd+1)];
        outlier_ub = [outliers(outlier_bd) outliers(end)];

        for iout = 1:length(outlier_lb)
            cleaned_inds(round((outlier_lb(iout):outlier_ub(iout))+win/2)) = 0;
        end
    end
    cleaned_inds = logical(cleaned_inds);
    cleaned_signal = mysignal;
    c_signal(~cleaned_inds,i) = nan;
 
end


%% Plot mean traces to compare mutation status
f8 = figure;

xfac = 1;
yfac = 1;
fontsize = 6;

setFigure(f8,xfac,yfac,fontsize)

celltype_str = {'MCF10A-WT','MCF10A-AKTspec','MCF10A-ERKspec','184A1-WT','184A1-AKTspec','184A1-ERKspec'};
ligand_name = {'IGF','HRG','HGF','EGF','BTC','NS'};
colmap = lines(3);
linewidth = 1;
legh = [];
legstr = {};
for isite = sites_all
    sprop = siteprop(isite,extension);
    icol = strmatch(sprop.celltype,celltype_str,'exact');
    irow = sprop.lig_index;
    
    if irow < 6
        subplot(6,4,(irow-1)*4+(icol>3)*2+(sprop.lig_dose==100)+1)
    else
        subplot(6,4,(irow-1)*4+(icol>3)*2+1)
    end
    
    mean_sig = nanmean(c_signal(:,celltype == isite),2);
    std_sig = nanstd(c_signal(:,celltype == isite),[],2);
    
    hold on
    if irow == 1 && icol < 4 && sprop.lig_dose==20
        legh = [legh plot(timestamp,mean_sig,'Color',colmap(mod(icol-1,3)+1,:))];
        legstr{end+1} = sprop.celltype(8:end);
    else
        plot(timestamp,mean_sig,'Color',colmap(mod(icol-1,3)+1,:))
    end
    
    tmpx = [timestamp; flipud(timestamp)];
    tmpy = [mean_sig + 2*std_sig/sqrt(sum(celltype == isite)); flipud(mean_sig - 2*std_sig/sqrt(sum(celltype == isite)))];
    
    ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp, 'FaceColor', colmap(mod(icol-1,3)+1,:)*0.1+0.9, 'EdgeColor', colmap(mod(icol-1,3)+1,:)*0.1+0.9,'LineWidth',linewidth);
    ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', colmap(mod(icol-1,3)+1,:)*0.3+0.7,'LineWidth',linewidth);
    
    plot(get(gca,'XLim'),[0 0],'--k')
    
    set(gca,'XLim',[50 800])
    
    plot([120 120],[-0.01 0.02],'b--')
    
    if irow == 6 && icol == 1
        xlabel('time [min]')
    end
    if irow == 1 && icol == 1 && sprop.lig_dose==20
        ylabel('log_{10} FOXO3a Cyt/Nuc ratio');
    end
    set(gca,'YLim',[-0.02 0.025])
    
    title(sprintf('%s %g - %s',sprop.lig_name,sprop.lig_dose,sprop.celltype(1:(5+(icol<4)))))
end
legend(legh,legstr)


%% Plot pulsing scores to compare mutation status
puls_thres = [.6 .35];

celltype_str = {'MCF10A-WT','MCF10A-AKTspec','MCF10A-ERKspec','184A1-WT','184A1-AKTspec','184A1-ERKspec'};
ligand_name = {'IGF','HRG','HGF','EGF','BTC','NS'};
celltype_two = {'MCF10A','184A1'};

ratio_puls_mat = nan(2*5+1,3,2);
ratio_puls_lb = nan(2*5+1,3,2);
ratio_puls_ub = nan(2*5+1,3,2);
ratio_fun1 = @(x) sum(x > puls_thres(1)) / length(x);
ratio_fun2 = @(x) sum(x > puls_thres(2)) / length(x);

legstr = cell(1,size(ratio_puls_mat,2));
xlab = cell(1,size(ratio_puls_mat,1));
for isite = sites_all
    sprop = siteprop(isite,extension);
    icol = 2*(strmatch(sprop.lig_name,ligand_name,'exact')-1)+(sprop.lig_dose==100)+1;
    irow = mod(strmatch(sprop.celltype,celltype_str,'exact')-1,3)+1;
    icell = (strmatch(sprop.celltype,celltype_str,'exact')>3)+1;
    
    ratio_puls_mat(icol,irow,icell) = nansum(dists(celltype == isite) > puls_thres(icell)) / sum(celltype == isite);
    if sum(celltype == isite) > 1
        if icell == 1
            tmp = bootci(2000,{ratio_fun1,dists(celltype == isite)},'alpha',.32);
        else
            tmp = bootci(2000,{ratio_fun2,dists(celltype == isite)},'alpha',.32);
        end
        ratio_puls_lb(icol,irow,icell) = tmp(1);
        ratio_puls_ub(icol,irow,icell) = tmp(2);
    end

    if irow == 1 && icell == 1
        xlab{icol} = sprintf('%s %g',sprop.lig_name,sprop.lig_dose);
    end
    
    if icol == 1 && icell == 1
        legstr{irow} = sprop.celltype(8:end);
    end
end

f8 = figure;

xfac = 1;
yfac = .6;
fontsize = 6;

setFigure(f8,xfac,yfac,fontsize)

for icell = 1:size(ratio_puls_mat,3)
    subplot(1,2,icell)
    
    h = bar(ratio_puls_mat(:,:,icell));
    set(h(1),'FaceColor',[0 0 1])
    set(h(2),'FaceColor',[1 0 0])
    set(h(3),'FaceColor',[1 0 1])
    x1 = get(get(h(1),'children'), 'xdata');
    x1 = mean(x1([1 3],:));
    x2 = get(get(h(2),'children'), 'xdata');
    x2 = mean(x2([1 3],:));
    x3 = get(get(h(3),'children'), 'xdata');
    x3 = mean(x3([1 3],:));
    hold on
    h=errorbar(x1,ratio_puls_mat(:,1,icell),ratio_puls_mat(:,1,icell)-ratio_puls_lb(:,1,icell),ratio_puls_ub(:,1,icell)-ratio_puls_mat(:,1,icell),'k'); set(h,'linestyle','none')
    h=errorbar(x2,ratio_puls_mat(:,2,icell),ratio_puls_mat(:,2,icell)-ratio_puls_lb(:,2,icell),ratio_puls_ub(:,2,icell)-ratio_puls_mat(:,2,icell),'k'); set(h,'linestyle','none')
    h=errorbar(x3,ratio_puls_mat(:,3,icell),ratio_puls_mat(:,3,icell)-ratio_puls_lb(:,3,icell),ratio_puls_ub(:,3,icell)-ratio_puls_mat(:,3,icell),'k'); set(h,'linestyle','none')
    
    set(gca,'XTick',1:size(ratio_puls_mat,1),'XTickLabel',xlab,'XLim',[.5 size(ratio_puls_mat,1)+.5])
    set(gca,'YLim',[0 .65])
    title(celltype_two{icell})
    ylabel('ratio pulsing cells')
end
legend(legstr)


%% Plot mean dynamics of traces with puls_score < .1 against those with puls_score > .1
f5 = figure;

xfac = 1;
yfac = 1;
fontsize = 4;

setFigure(f5,xfac,yfac,fontsize)

rowstocols = 0.5;
nrows = 6;
ncols = 12;

for isite = sites_all
    subplot(nrows,ncols,subplotpos(isite,ncols))
    hold on
    
    plot(timestamp,nanmean(c_signal(:,celltype==isite & dists < .1),2),'Color','r')
    plot(timestamp,nanmean(c_signal(:,celltype==isite & dists >= .1),2),'Color','k')
    
    set(gca,'XLim',[0 800],'YLim',[-.015 .025])
    
    sprop = siteprop(isite,extension);
    title(sprintf('%s %g',sprop.lig_name,sprop.lig_dose))
end

%% Plot single traces for specific condition
figure

mysite = 51;

nrows = 10;
ncols = 10;

myinds = find(celltype == mysite & dists >= .6);
% myinds = find(celltype == mysite & dists < .1);

for isig = 1:min([length(myinds) nrows*ncols])
    subplot(nrows,ncols,isig)
    plot(timestamp,c_signal(:,myinds(isig)))
    set(gca,'XLim',[0 800],'YLim',[-.04 .04])
end

%% Compare single traces for 2 specific conditions
mysites = [49 51];

nrows = 10;
ncols = 10;

[dists_sorted ind_sorted] = sort(dists,'descend');
nsubplot = min([sum(celltype == mysites(1)); sum(celltype == mysites(2)); nrows*ncols]);
c_signal_sorted = c_signal(:,ind_sorted);

for isite = mysites
    f = figure;

    xfac = 1;
    yfac = 1;
    fontsize = 4;

    setFigure(f,xfac,yfac,fontsize)
    
    myind = round(linspace(1,sum(celltype == isite),nsubplot));
    mysig = find(celltype(ind_sorted) == isite);
    
    for isig = 1:length(myind)
        subplot(nrows,ncols,isig)
        plot(timestamp,c_signal_sorted(:,mysig(myind(isig))))
        title(sprintf('PulsScore = %g',dists_sorted(mysig(myind(isig)))))
        set(gca,'XLim',[0 800],'YLim',[-.04 .04])
    end
end

%% Compare pulsing histograms accross ligands
% highdoses = [1 24 25 48 49 66]; % MCF10AWT
% highdoses = [3 22 27 46 51 65]; % MCF10AAKTspec
% highdoses = [5 20 29 44 53 64]; % MCF10AERKspec
% highdoses = [7 18 31 42 55 63]; % 184A1WT
highdoses = [9 16 33 40 57 62]; % 184A1AKTspec
% highdoses = [11 14 35 38 59 61]; % 184A1ERKspec

% kswidth = 15;
kswidth = 5;

highinds = ismember(celltype,highdoses);
dists_boxcox = boxcox(dists(highinds));

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
set(gca,'XLim',[-6 2])
title(s.celltype)
legend(legh,legstr)

xlabel('pulsatory score')
ylabel('probability density')