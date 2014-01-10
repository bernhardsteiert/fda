%% Generate table for Google Spreadsheet

close all
clear all
clc

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
addpath(grabdataPath)

% Get properties of sites by calling siteprop(site)
% All ligands + EGF with Inhibitors
% sites = [1 2 4:10 17:-1:11 24:30 37:-1:31 41 42 44:50 57:-1:51 64:69];
sites_for_harmonics = [4:10 17:-1:11 24:30 37:-1:31 44:50 57:-1:51 64:70];
sites = sites_for_harmonics;

% sites = [1 4];

times = cell(0);
signals = cell(0);
celltype = [];

for isite = sites
    if exist(remotepath,'dir')
        [times{end+1},intensity] = grabdata(isite);
    else
        load(['./Workspaces/site_' num2str(isite)])
        times{end+1} = timestamp;
    end

    log_trafo = 1; % log-transform signal

    if log_trafo
        signals{end+1} = log10(intensity);
    else
        signals{end+1} = intensity;
    end
    
    celltype = [celltype ones(1,size(intensity,2))*isite];
end

timestamp = times{1}; % same time sampling for all data sets
c_signal = cell2mat(signals);

close all

time_range = [50 200];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;


%% Generate spline fits to data-sets given in sites_for_harmonics
% close all

nbasis = 20;
% time_range = [min(timestamp) max(timestamp)];

ind_harm = ismember(celltype,sites_for_harmonics);
ind_fit = ~ind_harm;

basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
smoothed_data = smooth_basis(timestamp(range_ind),c_signal(range_ind,ind_harm),basis);


%% Make FPCA with data generated in previous block
% close all

nharm = 3;
% c_signal_pcastr = pca_fd(smoothed_data, nharm, fdPar(basis, int2Lfd(2), 0));
c_signal_pcastr = pca_fd(smoothed_data, nharm, fdPar(basis, int2Lfd(2), 0), 0); % WITHOUT CENTERING!!
% c_signal_pcastr = varmx_pca(c_signal_pcastr);


%% Generate table
close all

% Define principal components to be plotted
pcs = [2 3];

angle1 = -25;
angle2 = 5;
angle3 = -5;

Rmat1 = [cos(2*pi*angle1/360)   sin(2*pi*angle1/360)   0; ...
         -sin(2*pi*angle1/360)  cos(2*pi*angle1/360)   0; ...
         0                      0                      1];
     
Rmat2 = [cos(2*pi*angle2/360)   0  sin(2*pi*angle2/360); ...
         0                      1  0; ...
         -sin(2*pi*angle2/360)  0  cos(2*pi*angle2/360)];
     
Rmat3 = [1 0                      0; ...
         0 cos(2*pi*angle3/360)   sin(2*pi*angle3/360); ...
         0 -sin(2*pi*angle3/360)  cos(2*pi*angle3/360)];
     
Rmat = Rmat3 * Rmat2 * Rmat1;

flipharm = ones(1,nharm);
flipharm(1:3) = [1 1 -1]; % Unregistered

unitypes = unique(celltype(ind_harm));

site_lig_ind = [];
site_lig_name = {};
site_lig_dose = [];
site_inh_name = {};
site_inh_dose = [];

% Exclude Inhibitor data
for isite = sites_for_harmonics   % only plot ligands used for harmonics
    s = siteprop(isite);
    site_lig_ind = [site_lig_ind s.lig_index];
    site_lig_name{end+1} = s.lig_name;
    site_lig_dose = [site_lig_dose s.lig_dose];
    site_inh_name{end+1} = s.inh_name;
    site_inh_dose = [site_inh_dose s.inh_dose];
end

sites_remain = find(~site_inh_dose);
uni_lig = unique(site_lig_ind(sites_remain));
lig_min = min(site_lig_dose(sites_remain));
lig_max = max(site_lig_dose(sites_remain));


flipped_scores = repmat(flipharm,size(c_signal_pcastr.harmscr,1),1).*c_signal_pcastr.harmscr;

x_scores = (Rmat(pcs(1),:) * flipped_scores(:,1:3)')';
y_scores = (Rmat(pcs(2),:) * flipped_scores(:,1:3)')';

x_scores = (x_scores - min(x_scores))./range(x_scores); % Make range of x_scores = 1
y_scores = (y_scores - min(y_scores))*range(prctile(x_scores,[5 95]))./range(prctile(y_scores,[5 95]));

plot(x_scores,y_scores,'.','Color',[.7 .7 .7]);
axis equal

unstim = find(site_lig_dose == 0);
unstim = unstim(1:end-1); % BTC just a copy of FGF
x_scores_unstim = [];
y_scores_unstim = [];

for iunstim = unstim
    x_scores_unstim = [x_scores_unstim; x_scores(celltype(ind_harm) == sites_for_harmonics(iunstim))];
    y_scores_unstim = [y_scores_unstim; y_scores(celltype(ind_harm) == sites_for_harmonics(iunstim))];
end

x_med_unstim = median(x_scores_unstim);
y_med_unstim = median(y_scores_unstim);

xls_cell = {'Ligand','Concentration','PC2 center','PC3 center','Circle Size'};
filename = [datestr(now, 30) '_fPCA_bubbles_for_google_spreadsheet.csv'];

for isite = sites_remain
    ligname = site_lig_name{isite}; % Col 1
    concentration = site_lig_dose(isite); % Col 2
    
    x_med = median(x_scores(celltype(ind_harm) == sites_for_harmonics(isite),:)); % Col 3
    y_med = median(y_scores(celltype(ind_harm) == sites_for_harmonics(isite),:)); % Col 4
    
    rad_dist_sorted = sort(sqrt((x_scores(celltype(ind_harm) == sites_for_harmonics(isite),:) - x_med).^2 + (y_scores(celltype(ind_harm) == sites_for_harmonics(isite),:) - y_med).^2));
    circle_size = rad_dist_sorted(round(length(rad_dist_sorted)*.68268)); % Col 5
    
%     Strange: Radius and Angle have to be resized again for circle size, like PC2/3 ...
%     radius = median(sqrt((x_scores(celltype(ind_harm) == sites_for_harmonics(isite),:) - x_med_unstim).^2 + (y_scores(celltype(ind_harm) == sites_for_harmonics(isite),:) - y_med_unstim).^2)); % Col 6
%     resp_angle = atan((x_scores(celltype(ind_harm) == sites_for_harmonics(isite),:) - x_med_unstim)./(y_scores(celltype(ind_harm) == sites_for_harmonics(isite),:) - y_med_unstim));
%     resp_angle(resp_angle < 0) = pi + resp_angle(resp_angle < 0);
%     angle = median(resp_angle); % Col 7

    xls_cell(end+1,:) = {ligname(1:3),concentration,x_med,y_med,circle_size};
end

%% Write CSV File
xls_cell = xls_cell';

f = fopen(filename,'a+');
for icell = 1:length(xls_cell(:))
    if mod(icell-1,5) == 0 || icell <= 5
        fprintf(f,'%s',xls_cell{icell});
    else
        fprintf(f,'%f',xls_cell{icell});
    end
    
    if mod(icell,5) == 0
        fprintf(f,'\n');
    else
        fprintf(f,';');
    end
    
end
fclose(f);


return
% set(gca,'XLim',[min(x_scores) max(x_scores)]*1.1)
set(gca,'XLim',[min(x_scores) max(x_scores)]-.02)
aspRatioFig = posFig(3)/posFig(4);
posSubplot = get(gca,'Position');
aspRatioSubplot = posSubplot(3)/posSubplot(4);
% set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)) % PC1 vs PC3
%     set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)+.03) % PC1 vs PC2
set(gca,'YLim',[-.06 .12])

% if ilig_plot == 1
    ylabel(['PC ' num2str(pcs(2))])
%     arrow([-.13 -.1],[.2 -.1],'Width',.5,'Length',7)
%     set(gca,'YTick',-.1:.1:.2)
% end
% if ilig_plot == length(ligs_to_plot)-1
    xlabel(['PC ' num2str(pcs(1))])
%     set(gca,'YTickLabel',[])
% end
% if ilig_plot == length(ligs_to_plot)
%     arrow([-.13 .08],[.05 .2],'Width',.5,'Length',7)
%     set(gca,'YTickLabel',[])
% end

% lline = .13;
% plot([x_med x_med],get(gca,'YLim')*.95,'k-','LineWidth',2)
% plot(x_med,y_med,'kx','LineWidth',3)
% alpha = pi/4;
% angles = linspace(0,alpha,51);
% plot([x_med x_med+lline*sin(pi/4)],[y_med lline*cos(pi/4)],'k-','LineWidth',2)
% plot(x_med+lline*sin(angles)/2,y_med+lline*cos(angles)/2,'k-','LineWidth',2)


% set(gca,'children',flipud(get(gca,'children')))

h = get(gca);
subplot(nrows,ncols,[4 10])
clim = [0 1];
set(gca,'CLim',clim)
colormap(colmap)
% colorbar('Location','North','XTick',log10([1 2.5 5 10 20 50 100]),'XTickLabel',[0 2.5 5 10 20 50 100]) % Horizontal colorbar
colorbar('Location','West','YTick',linspace(range(clim)./(2*length(ligs_to_plot)),1-range(clim)./(2*length(ligs_to_plot)),length(ligs_to_plot)),'YTickLabel',leg_str, 'TickLength', [0 0]) % Vertical colorbar
set(gca,'Visible','off')
% text(0,-.05,'Ligand')
possubp = get(gca,'Position');
set(gca,'Position',[possubp(1)-.03 possubp(2:4)])

resp_strength = [];
resp_angle = [];
for isite = ligs_to_plot(end:-1:1)
    tmpind = find(site_lig_ind(sites_remain) == uni_lig(isite));
    iplot = find(site_lig_dose(tmpind) == 100); % Only high dose
    
    resp_strength = padconcatenation(resp_strength,sqrt((x_scores(celltype(ind_harm) == sites_for_harmonics(tmpind(iplot)),:) - x_med).^2 + (y_scores(celltype(ind_harm) == sites_for_harmonics(tmpind(iplot)),:) - y_med).^2),2);
    resp_angle = padconcatenation(resp_angle,atan((x_scores(celltype(ind_harm) == sites_for_harmonics(tmpind(iplot)),:) - x_med)./(y_scores(celltype(ind_harm) == sites_for_harmonics(tmpind(iplot)),:) - y_med)),2);
end
resp_angle(resp_angle < 0) = pi + resp_angle(resp_angle < 0);

% Response strength plot
subplot(nrows,ncols,5)
boxplot(resp_strength)
possubp1 = get(gca,'Position');
% title('Response strength')
xlabel('Response strength')
hold on
set(gca,'XTick',1:6,'XTickLabel',leg_str)
midpoint = .1;
range_bars = .07;
plot(get(gca,'XLim'),[1 1]*midpoint,'k--')
plot(get(gca,'XLim'),[1 1]*(midpoint+range_bars),'k:')
plot(get(gca,'XLim'),[1 1]*(midpoint-range_bars),'k:')
% set(gca,'YLim',[pi/6 2/3*pi]) % Almost everything visible
set(gca,'YLim',[0 .2],'YTick',[midpoint-range_bars midpoint midpoint+range_bars],'YTickLabel',{'none','medium','high'})
set(gca,'Position',[.63 possubp1(2) .14 possubp1(4)])

% Response strength plot
sps = subplot(nrows,ncols,6);
hold on
legendstyles = [];
for ip = 1:size(resp_strength,2)

    baredges = linspace(0,.2,6);
    barheight = histc(resp_strength(:,ip),baredges)./sum(~isnan(resp_strength(:,ip)));
    legendstyles = [legendstyles plot(baredges,barheight,'Color',colmap(ip,:),'LineWidth',2)];
%     [f,xi] = ksdensity(resp_strength(:,ip));
%     legendstyles = [legendstyles plot(xi,f,'Color',colmap(ip,:))];

end
% title('Response strength')
ylabel('distribution')

% legend(sps,legendstyles,leg_str)

% set(gca,'XTick',1:6,'XTickLabel',leg_str)
midpoint = .1;
range_bars = .07;
% plot([1 1]*midpoint,get(gca,'YLim'),'k--')
% plot([1 1]*(midpoint+range_bars),get(gca,'YLim'),'k:')
% plot([1 1]*(midpoint-range_bars),get(gca,'YLim'),'k:')
% % set(gca,'YLim',[pi/6 2/3*pi]) % Almost everything visible
set(gca,'XLim',[0 .2],'XTick',[midpoint-range_bars midpoint midpoint+range_bars],'XTickLabel',{'none','medium','high'})
possubp1 = get(gca,'Position');
set(gca,'Position',[possubp1(1)+.04 possubp1(2) possubp1(3)+.03 possubp1(4)])

% Response angle plot
subplot(nrows,ncols,11)
boxplot(resp_angle)
possubp2 = get(gca,'Position');
xlabel('Response angle')
hold on
set(gca,'XTick',1:6,'XTickLabel',leg_str)
% 3*pi/8 is middle??
% midpoint = 7*pi/16;
midpoint = 1.3;
range_bars = 0.3;
plot(get(gca,'XLim'),[1 1]*midpoint,'k--')
plot(get(gca,'XLim'),[1 1]*(midpoint+range_bars),'k:')
plot(get(gca,'XLim'),[1 1]*(midpoint-range_bars),'k:')
% set(gca,'YLim',[pi/6 2/3*pi]) % Almost everything visible
set(gca,'YLim',[.8 1.8],'YTick',[midpoint-range_bars midpoint midpoint+range_bars],'YTickLabel',{'transient','both','persistent'})
% possubp2 = get(gca,'Position');
set(gca,'Position',[.63 possubp2(2) .14 possubp2(4)])
% set(gca,'Position',possubp2-.01)

% Response angle plot
sps2 = subplot(nrows,ncols,12);
% boxplot(resp_angle)
hold on
legendstyles = [];
for ip = 1:size(resp_angle,2)

    baredges = linspace(0.8,1.8,6);
    barheight = histc(resp_angle(:,ip),baredges)./sum(~isnan(resp_angle(:,ip)));
    legendstyles = [legendstyles plot(baredges,barheight,'Color',colmap(ip,:),'LineWidth',2)];
%     [f,xi] = ksdensity(resp_strength(:,ip));
%     legendstyles = [legendstyles plot(xi,f,'Color',colmap(ip,:))];

end

% title('Response angle')
ylabel('distribution')
% hold on
% set(gca,'XTick',1:6,'XTickLabel',leg_str)
% % 3*pi/8 is middle??
% % midpoint = 7*pi/16;
midpoint = 1.3;
range_bars = 0.3;
% plot(get(gca,'XLim'),[1 1]*midpoint,'k--')
% plot(get(gca,'XLim'),[1 1]*(midpoint+range_bars),'k:')
% plot(get(gca,'XLim'),[1 1]*(midpoint-range_bars),'k:')
% % set(gca,'YLim',[pi/6 2/3*pi]) % Almost everything visible
set(gca,'XLim',[.8 1.8],'XTick',[midpoint-range_bars midpoint midpoint+range_bars],'XTickLabel',{'transient','both','persistent'})
possubp2 = get(gca,'Position');
set(gca,'Position',[possubp2(1)+.04 possubp2(2) possubp2(3)+.03 possubp2(4)])

return

%% Plot EGF early vs. late time-points
close all

time_range = [50 510];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

f = figure;

posFig = get(gcf,'Position');
posFig(4) = posFig(4)/1.5;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./25);

nrows = 3;
ncols = 5;
subplot(nrows,ncols,[1 2 3 6 7 8 11 12 13])

site_lig_ind = [];
site_lig_name = {};
site_lig_dose = [];
site_inh_name = {};
site_inh_dose = [];

% Exclude Inhibitor data
for isite = sites_for_harmonics   % only plot ligands used for harmonics
    s = siteprop(isite);
    site_lig_ind = [site_lig_ind s.lig_index];
    site_lig_name{end+1} = s.lig_name;
    site_lig_dose = [site_lig_dose s.lig_dose];
    site_inh_name{end+1} = s.inh_name;
    site_inh_dose = [site_inh_dose s.inh_dose];
end

iplot = 1; % EGF high dose
first_n = 15; % Plot first_n traces colored
set(f,'DefaultAxesColorOrder',lines(first_n))

plotsignals = find(celltype == sites_for_harmonics(iplot));
plot(repmat(timestamp(range_ind),1,length(plotsignals)),c_signal(range_ind,plotsignals),'g','color',[0.7 0.7 0.7])
hold on
plot(repmat(timestamp(range_ind),1,first_n),c_signal(range_ind,plotsignals(1:first_n)))
plot(timestamp(range_ind),nanmean(c_signal(range_ind,plotsignals),2),'color','k','LineWidth',2)
title([site_lig_name{iplot} num2str(site_lig_dose(iplot)) ' [ng/ml]'])

xlabel('time [min]')
ylabel('log_{10} FOXO3a Cyt/Nuc ratio');

ylim = [-1 1]*.04;
if ~log_trafo
    ylim = 10.^ylim;
end
plot([200 200],ylim,'k--')
text(75,.035,'early response\newline(deterministic)')
text(300,.035,'late response\newline(stochastic)')

set(gca,'XLim',time_range,'YLim',ylim)
set(gca,'XTick',50:50:500)

% Figure 2B: Eigenfunctions (new - rotated)
pcs = 1:3;

% rowstocols = 1;
% nrows = ceil(nharm^rowstocols);
% ncols = ceil(nharm / nrows);

time_range = [50 200];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;
times_fine = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),501);

harm_eval = repmat(flipharm(1:nharm),length(times_fine),1) .* eval_fd(c_signal_pcastr.harmfd,times_fine);
harm_eval_rescale = 2 * repmat(sqrt(c_signal_pcastr.values(1:nharm))',length(times_fine),1) .* harm_eval;
% <-- Why does Ramsay do that?? Destroys normalization and may not be done if harmonics are rotated!
harm_eval_rescale = harm_eval;

% posFig = get(gcf,'Position');
% posFig(3) = posFig(3)/2.5;
% set(gcf,'Position',posFig)
% set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./20);

harmscr = (Rmat * flipped_scores(:,1:3)')';
propvar = squeeze(sum(harmscr.^2));
propvar = propvar./sum(propvar);
propvar = propvar.*sum(c_signal_pcastr.varprop(1:3));

for iplot = 1:nharm
    subplot(nrows,ncols,ncols*(iplot-1)+ncols-1)
    
    if iplot <= length(pcs)
        tmpplot = sum(repmat(Rmat(iplot,:),size(harm_eval_rescale,1),1) .* harm_eval_rescale(:,pcs),2);
        varper = propvar(iplot);
    else
        tmpplot = harm_eval_rescale(:,iplot);
        varper = c_signal_pcastr.varprop(iplot);
    end
    
    plot(times_fine,tmpplot)
    title(['Harmonic ' num2str(iplot)])
    set(gca,'XLim',time_range)
    set(gca,'YLim',[min(min(harm_eval_rescale)) max(max(harm_eval_rescale))]*1.2)
    
    hold on
    plot(time_range,[0 0],'--')
    
%     title(['Variance explained: ' num2str(varper*100,3) '%'])
    
    if iplot == nharm
%         xlabel('time [min]')
        set(gca,'XTick',50:50:200)
    else
        set(gca,'XTick',[])
    end
    
    set(gca,'YTick',[])
end

load('harm_basis.mat') % Contains only harm_basis from all data-sets


time_range = [200 510];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;
times_fine = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),501);

basis_eval = eval_basis(harm_basis,times_fine);

for iplot = 1:size(basis_eval,2)
    subplot(nrows,ncols,ncols*(iplot-1)+ncols)
    
    tmpplot = basis_eval(:,iplot);
    
    plot(times_fine,tmpplot)
    title(['Harmonic ' num2str(iplot)])
    set(gca,'XLim',time_range)
    set(gca,'YLim',[min(min(basis_eval)) max(max(basis_eval))]*1.2)
    
    hold on
    plot(time_range,[0 0],'--')
    
%     title(['Variance explained: ' num2str(varper*100,3) '%'])
    
    if iplot == nharm
%         xlabel('time [min]')
        set(gca,'XTick',200:100:500)
    else
        set(gca,'XTick',[])
    end
    
    set(gca,'YTick',[])
end