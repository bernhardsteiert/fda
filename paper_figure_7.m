%% Paper argumentation for fPCA part (figure 2)
% 0. Cut signals in 2 parts dominated by a) Deterministic rise (t = [50 200]) b) Stochastic oscillations (t = [200 650])
% 1. Do fPCA for all ligands at all doses (no inhibitors) for t = [50 200]
%    --> PC2vs3 - IGF one axis - BTC other axis - EGF/HGF/HRG inbetween - EPR/FGF not responding
%    Plot first 4 harmonics --> argument for choosing PCs 2 and 3
% 2. Visualize this with ellipses instead of points for IGF and BTC

% The following points are moved to figure 3
% 3. Subtract mean + first 3 harmonics (plot) from data for t = [200 650]
%    --> Histograms for no dose + highest dose with threshold
%    cell-to-cell heterogeneity <-> time-courses to show oscillations
%    Ligands triggering subpopulation to oscillate; Others don't
% 4. Do same analysis for inhibitor data (EGF and HGF): MEKi induces oscillations; AKTi represses oscillations
%    (Counter-intuitive: EGF-MEKi population mean like IGF, oscillations however not suppressed as for IGF)


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

% return


%% Register to mean value in time range
close all

time_range = [50 200];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

time_range_reg = [50 120];

[tmp range_ind_min] = min(abs(timestamp - time_range_reg(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range_reg(2)));
range_ind_reg = range_ind_min:range_ind_max;

figure

posFig = get(gcf,'Position');
% posFig(4) = posFig(4)/2;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./10);

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

nrows = 7;
ncols = 7;

for iplot = 1:length(sites_for_harmonics)
    subplot(nrows,ncols,iplot)
    
%     c_signal(:,celltype == sites_for_harmonics(iplot)) = c_signal(:,celltype == sites_for_harmonics(iplot)) - nanmean(nanmean(c_signal(range_ind_reg,celltype == sites_for_harmonics(iplot))));
    
    plot(repmat(timestamp(range_ind),1,sum(celltype == sites_for_harmonics(iplot))),c_signal(range_ind,celltype == sites_for_harmonics(iplot)),'g','color',[0.7 0.7 0.7])
    hold on
    plot(timestamp(range_ind),nanmean(c_signal(range_ind,celltype == sites_for_harmonics(iplot)),2),'color','k','LineWidth',2)
    title([site_lig_name{iplot} num2str(site_lig_dose(iplot))])
    
    ylim = [-1 1]*.04;
    if ~log_trafo
        ylim = 10.^ylim;
    end
    set(gca,'XLim',time_range,'YLim',ylim)

end

% c_signal = c_signal - repmat(nanmean(c_signal(range_ind_reg,:),1),size(c_signal,1),1);


%% Generate spline fits to data-sets given in sites_for_harmonics
% close all

nbasis = 20;
% time_range = [min(timestamp) max(timestamp)];

ind_harm = ismember(celltype,sites_for_harmonics);
ind_fit = ~ind_harm;

basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
smoothed_data = smooth_basis(timestamp(range_ind),c_signal(range_ind,ind_harm),basis);

% return

% f = figure;
% set(f,'DefaultAxesColorOrder',jet(size(c_signal(1,ind_harm),2)))
% hold on
% 
% plot(smoothed_data)
% plot(timestamp(range_ind),c_signal(range_ind,ind_harm),'o')

%% Make FPCA with data generated in previous block
% close all

nharm = 3;
% c_signal_pcastr = pca_fd(smoothed_data, nharm, fdPar(basis, int2Lfd(2), 0));
c_signal_pcastr = pca_fd(smoothed_data, nharm, fdPar(basis, int2Lfd(2), 0), 0); % WITHOUT CENTERING!!
% c_signal_pcastr = varmx_pca(c_signal_pcastr);

% return

% plot_pca_fd(c_signal_pcastr, 1, 0)

% c_signal_rotpcastr = varmx_pca(c_signal_pcastr);
% plot_pca_fd(c_signal_rotpcastr, 1, 0)

%% Figure 2A: fPCA with dose-dependend colored ligand level
close all

% Define principal components to be plotted
pcs = [2 3];

% Unregistered data - Harm 1: DC; Harm 2: Rise; Harm 3: Peak
% angle1 = 20;
% angle2 = 0;
% angle3 = -10;
angle1 = -25;
angle2 = 5;
angle3 = -5;
% Registered data - Harm 1: SS_afterStim; Harm 2: SS_preStim; Harm 3: Peak
% angle1 = 15;
% angle2 = -7;
% Old:
% angle = -10; % rotation angle to right [degree]

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
% Has to be adjusted according to angles
% flipharm(1:4) = [1 1 -1 1]; % Unregistered
flipharm(1:3) = [1 1 -1]; % Unregistered

unitypes = unique(celltype(ind_harm));

% ligs_to_plot = 1:length(uni_lig); % Plot all
ligs_to_plot = [2 4 5 1 7 6]; % Plot all but FGF
% ligs_to_plot = [2 1 7];

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
ncolor = length(ligs_to_plot)-1;
colmap = flipud(hsv(ncolor+1));
% colmap = flipud(winter(ncolor));
color_doses = 10.^linspace(max(log10([lig_min 1])),log10(lig_max),ncolor);

rowstocols = 0.3;
% rowstocols = 0;
nrows = 2;
ncols = 5;

figure

posFig = get(gcf,'Position');
posFig(4) = posFig(4)/1.5;
% posFig(4) = posFig(4)/3;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./25);

flipped_scores = repmat(flipharm,size(c_signal_pcastr.harmscr,1),1).*c_signal_pcastr.harmscr;

STD = 1;                     %# 2 standard deviations
conf = 2*normcdf(STD)-1;     %# covers around 95% of population
scale = chi2inv(conf,2);     %# inverse chi-squared with dof=#dimensions
    
subplot(nrows,ncols,[1 2 3 6 7 8])
box on
hold on

x_scores = (Rmat(pcs(1),:) * flipped_scores(:,1:3)')';
y_scores = (Rmat(pcs(2),:) * flipped_scores(:,1:3)')';
plot(x_scores,y_scores,'.','Color',[.7 .7 .7]);

unstim = find(site_lig_dose == 0);
unstim = unstim(1:end-1); % BTC just a copy of FGF
x_scores_unstim = [];
y_scores_unstim = [];

for iunstim = unstim
    x_scores_unstim = [x_scores_unstim; x_scores(celltype(ind_harm) == sites_for_harmonics(iunstim))];
    y_scores_unstim = [y_scores_unstim; y_scores(celltype(ind_harm) == sites_for_harmonics(iunstim))];
end

x_med = median(x_scores_unstim);
y_med = median(y_scores_unstim);

leg_str = {};
color_ind = 1;

for isite = ligs_to_plot(end:-1:1)
    tmpind = find(site_lig_ind(sites_remain) == uni_lig(isite));
    leg_str{end+1} = site_lig_name{sites_remain(tmpind(1))};
    iplot = find(site_lig_dose(tmpind) == 100); % Only high dose
    
    % Colored
    mycolor = colmap(color_ind,:);
    color_ind = color_ind + 1;
    plot(x_scores(celltype(ind_harm) == sites_for_harmonics(tmpind(iplot)),:),y_scores(celltype(ind_harm) == sites_for_harmonics(tmpind(iplot)),:),'o','MarkerFaceColor',mycolor,'MarkerEdgeColor','none','MarkerSize',4);

    %# substract mean
    Mu = mean( [x_scores(celltype(ind_harm) == sites_for_harmonics(tmpind(iplot)),:) y_scores(celltype(ind_harm) == sites_for_harmonics(tmpind(iplot)),:)] );
    X0 = bsxfun(@minus, [x_scores(celltype(ind_harm) == sites_for_harmonics(tmpind(iplot)),:) y_scores(celltype(ind_harm) == sites_for_harmonics(tmpind(iplot)),:)], Mu);

    %# eigen decomposition [sorted by eigen values]
    Cov = cov(X0) * scale;
    [V D] = eig(Cov);
    [D order] = sort(diag(D), 'descend');
    D = diag(D);
    V = V(:, order);

    t = linspace(0,2*pi,100);
    e = [cos(t) ; sin(t)];        %# unit circle
    VV = V*sqrt(D);               %# scale eigenvectors
    e = bsxfun(@plus, VV*e, Mu'); %#' project circle back to orig space

    %# plot cov and major/minor axes
%         plot(e(1,:), e(2,:), 'Color',mycolor);

    tmpx = e(1,:);
    tmpy = e(2,:);
%     ltmp = patch(tmpx, tmpy, ones(size(tmpx)), ones(size(tmpx)));
%     set(ltmp, 'FaceColor', mycolor, 'EdgeColor', 'none', 'FaceAlpha', 1);


end

set(gca,'XLim',[min(x_scores) max(x_scores)]*1.1)
aspRatioFig = posFig(3)/posFig(4);
posSubplot = get(gca,'Position');
aspRatioSubplot = posSubplot(3)/posSubplot(4);
set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)) % PC1 vs PC3
%     set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)+.03) % PC1 vs PC2

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

lline = .13;
plot([x_med x_med],get(gca,'YLim')*.95,'k-','LineWidth',2)
plot(x_med,y_med,'kx','LineWidth',3)
alpha = pi/4;
angles = linspace(0,alpha,51);
plot([x_med x_med+lline*sin(pi/4)],[y_med lline*cos(pi/4)],'k-','LineWidth',2)
plot(x_med+lline*sin(angles)/2,y_med+lline*cos(angles)/2,'k-','LineWidth',2)
% set(gca,'children',flipud(get(gca,'children')))

h = get(gca);
subplot(nrows,ncols,[4 9])
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
title('Response strength')
hold on
set(gca,'XTick',1:6,'XTickLabel',leg_str)
midpoint = .1;
range_bars = .07;
plot(get(gca,'XLim'),[1 1]*midpoint,'k--')
plot(get(gca,'XLim'),[1 1]*(midpoint+range_bars),'k:')
plot(get(gca,'XLim'),[1 1]*(midpoint-range_bars),'k:')
% set(gca,'YLim',[pi/6 2/3*pi]) % Almost everything visible
set(gca,'YLim',[0 .2],'YTick',[midpoint-range_bars midpoint midpoint+range_bars],'YTickLabel',{'none','medium','high'})
possubp1 = get(gca,'Position');
set(gca,'Position',[.75 possubp1(2) .24 possubp1(4)])

% Response angle plot
subplot(nrows,ncols,10)
boxplot(resp_angle)
title('Response angle')
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
possubp2 = get(gca,'Position');
set(gca,'Position',[.75 possubp2(2) .24 possubp2(4)])

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