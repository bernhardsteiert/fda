%% Paper argumentation for fPCA part (figure 2)
% 0. Cut signals in 2 parts dominated by a) Deterministic rise (t = [50 200]) b) Stochastic oscillations (t = [200 650])
% 1. Do fPCA for all ligands at all doses (no inhibitors) for t = [50 200]
%    --> PC2vs3 - IGF one axis - BTC other axis - EGF/HGF/HRG inbetween - EPR/FGF not responding
%    Plot first 4 harmonics --> argument for choosing PCs 2 and 3
% 2. Add EGF inhibitor data and see where it goes (probably: MEKi --> IGF; AKTi --> not responding)
%    OR: Do the same with different inhibitor concentrations
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
sites = [1 2 4:10 17:-1:11 24:30 37:-1:31 41 42 44:50 57:-1:51 64:69];
sites_for_harmonics = [4:10 17:-1:11 24:30 37:-1:31 44:50 57:-1:51 64:69];

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

return

%% Plot raw data of defined ligands
figure

plot_ligs = [1 4 17 31];
plot_name = {'EGF-MEKi','EGF','IGF','No Stim'};
nrows = 2;
ncols = 2;

for iplot = 1:length(plot_ligs)
    subplot(nrows,ncols,iplot)
    
    plot(repmat(timestamp(range_ind),1,sum(celltype == plot_ligs(iplot))),c_signal(range_ind,celltype == plot_ligs(iplot)),'g','color',[0.7 0.7 0.7])
    hold on
    plot(timestamp(range_ind),nanmean(c_signal(range_ind,celltype == plot_ligs(iplot)),2),'color','k','LineWidth',2)
    title(plot_name{iplot})
    
    set(gca,'YLim',[-1 1]*.04)

end


%% Generate spline fits to data-sets given in sites_for_harmonics
% close all

nbasis = 20;
% time_range = [min(timestamp) max(timestamp)];
time_range = [50 200];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

ind_harm = ismember(celltype,sites_for_harmonics);
ind_fit = ~ind_harm;

basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
smoothed_data = smooth_basis(timestamp(range_ind),c_signal(range_ind,ind_harm),basis);

return

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal(1,ind_harm),2)))
hold on

plot(smoothed_data)
plot(timestamp(range_ind),c_signal(range_ind,ind_harm),'o')

%% Make FPCA with data generated in previous block
% close all

nharm = 4;
% c_signal_pcastr = pca_fd(smoothed_data, nharm);
c_signal_pcastr = pca_fd(smoothed_data, nharm, fdPar(basis, int2Lfd(2), 0), 0); % WITHOUT CENTERING!!
% c_signal_pcastr = varmx_pca(c_signal_pcastr);

return

plot_pca_fd(c_signal_pcastr, 1, 0)

% c_signal_rotpcastr = varmx_pca(c_signal_pcastr);
% plot_pca_fd(c_signal_rotpcastr, 1, 0)

%% Figure 2A: fPCA with dose-dependend colored ligand level
close all

% Define principal components to be plotted
pcs = [2 3];

% angle = 0;
angle = -10; % rotation angle to right [degree]

Rmat = [1                     0                     0; ...
        0                     cos(2*pi*angle/360)  sin(2*pi*angle/360); ...
        0                     -sin(2*pi*angle/360) cos(2*pi*angle/360)];

flipharm = ones(1,nharm);
flipharm(1:4) = [-1 1 -1 1];

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
ncolor = 201;
colmap = flipud(jet(ncolor));
color_doses = 10.^linspace(max(log10([lig_min 1])),log10(lig_max),ncolor);

rowstocols = 0.3;
nrows = ceil((length(uni_lig)+1)^rowstocols);
ncols = ceil((length(uni_lig)+1) / nrows);

figure

posFig = get(gcf,'Position');
posFig(4) = posFig(4)/2;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./15);

flipped_scores = repmat(flipharm,size(c_signal_pcastr.harmscr,1),1).*c_signal_pcastr.harmscr;

for ilig = 1:length(uni_lig)
    
    subplot(nrows,ncols,ilig)
    hold on
    tmpind = find(site_lig_ind(sites_remain) == uni_lig(ilig));
    title(site_lig_name{sites_remain(tmpind(1))})
    
    if ilig == 1
        ylabel(['PC ' num2str(pcs(2))])
    end
    if ilig == length(uni_lig)
        xlabel(['PC ' num2str(pcs(1))])
    end
    
    x_scores = (Rmat(pcs(1),:) * flipped_scores(:,1:3)')';
    y_scores = (Rmat(pcs(2),:) * flipped_scores(:,1:3)')';
    plot(x_scores,y_scores,'.','Color',[.7 .7 .7]);
    
    for isite = tmpind
        % Colored
        [tmp color_ind] = min(abs(color_doses-site_lig_dose(sites_remain(isite))));
        mycolor = colmap(color_ind,:);
        plot(x_scores(celltype(ind_harm) == sites_for_harmonics(sites_remain(isite)),:),y_scores(celltype(ind_harm) == sites_for_harmonics(sites_remain(isite)),:),'.','Color',mycolor);
    end
    
    set(gca,'XLim',[min(x_scores) max(x_scores)]*1.1)
    aspRatioFig = posFig(3)/posFig(4);
    posSubplot = get(gca,'Position');
    aspRatioSubplot = posSubplot(3)/posSubplot(4);
%     set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)-.03) % PC1 vs PC3
%     set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)+.03) % PC1 vs PC2
end

h = get(gca);
subplot(nrows,ncols,ilig+1)
clim = log10([lig_min lig_max]);
clim(1) = max([clim(1) 0]);
set(gca,'CLim',clim)
colormap(flipud(jet(ncolor)))
colorbar('Location','North','XTick',log10([1 2.5 5 10 20 50 100]),'XTickLabel',[0 2.5 5 10 20 50 100])
set(gca,'Visible','off')

%% Figure 2B: Eigenfunctions (new - rotated)
close all
pcs = 1:3;

rowstocols = 1;
nrows = ceil(nharm^rowstocols);
ncols = ceil(nharm / nrows);

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;
times_fine = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),501);

harm_eval = repmat(flipharm(1:nharm),length(times_fine),1) .* eval_fd(c_signal_pcastr.harmfd,times_fine);
harm_eval_rescale = 2 * repmat(sqrt(c_signal_pcastr.values(1:nharm))',length(times_fine),1) .* harm_eval;

figure

for iplot = 1:nharm
    subplot(nrows,ncols,iplot)
    
    if iplot <= length(pcs)
        tmpplot = sum(repmat(Rmat(iplot,:),size(harm_eval_rescale,1),1) .* harm_eval_rescale(:,pcs),2);
    else
        tmpplot = harm_eval_rescale(:,iplot);
    end

    
    plot(times_fine,tmpplot)
    xlabel(['Harmonic ' num2str(iplot)])
    set(gca,'XLim',time_range)
    set(gca,'YLim',[min(min(harm_eval_rescale)) max(max(harm_eval_rescale))]*1.2)
    
    hold on
    plot(time_range,[0 0],'--')
end

%% Fit additional data with basis from fPCA
% close all

harm_basis = create_fd_basis(c_signal_pcastr.harmfd);
% mean_fit = eval_fd(c_signal_pcastr.meanfd,timestamp(range_ind));
% smoothed_additional = smooth_basis(timestamp(range_ind),c_signal(range_ind,ind_fit)-repmat(mean_fit,1,sum(ind_fit)),harm_basis);
smoothed_additional = smooth_basis(timestamp(range_ind),c_signal(range_ind,ind_fit),harm_basis);
fitcoef = getcoef(smoothed_additional);
fitcoef = repmat(flipharm,size(fitcoef,2),1)'.*fitcoef;

return

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal(1,ind_fit),2)))
hold on


plot(smoothed_additional+c_signal_pcastr.meanfd)
plot(timestamp(range_ind),c_signal(range_ind,ind_fit),'o')

% %% Check if data is fitted by harmonics
% close all
% 
% nperplot = 8;
% nsubplots = ceil(size(c_signal(1,ind_fit),2)./nperplot);
% rowstocols = 0.5;
% time_range = [50 200];
% 
% [tmp range_ind_min] = min(abs(timestamp - time_range(1)));
% [tmp range_ind_max] = min(abs(timestamp - time_range(2)));
% range_ind = range_ind_min:range_ind_max;
% times_fine = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),501);
% harm_fine = eval_basis(harm_basis,times_fine);
% 
% nrows = ceil(nsubplots^rowstocols);
% ncols = ceil(nsubplots / nrows);
% 
% fitcoef = getcoef(smoothed_additional);
% 
% fitcoef_rot = fitcoef;
% fitcoef_rot(pcs(1),:) = (Rmat(pcs(1),:) * fitcoef(1:3,:));
% fitcoef_rot(pcs(2),:) = (Rmat(pcs(2),:) * fitcoef(1:3,:));
% 
% harm_fine_rot = harm_fine;
% harm_fine_rot(:,1:3) = (Rmat*harm_fine(:,1:3)')';
% 
% data_fpca_repr = fitcoef_rot'*harm_fine_rot';
% mean_fine = eval_fd(c_signal_pcastr.meanfd,times_fine);
% 
% ind_fit_no = find(ind_fit);
% 
% for iplot = 1:nsubplots
%     subplot(nrows,ncols,iplot)
%     
%     inds = ((iplot-1)*nperplot+1):min([size(c_signal(1,ind_fit),2) iplot*nperplot]);
%     plot(times_fine,repmat(mean_fine,1,length(inds))+data_fpca_repr(inds,:)')
%     hold on
%     plot(timestamp(range_ind),c_signal(range_ind,ind_fit_no(inds)),'o')
% 
% end

%% Figure 2C: Look what happens with inhibitor data in PCA space
close all

pcs = [2 3];

ind_included_panel1 = [31 4 17 64]; % NoStim; EGF; IGF; BTC (high doses)
ind_included_panel2 = [4 1 2];   % EGF high dose; + MEKi; + AKTi

mean_fine = eval_fd(c_signal_pcastr.meanfd,times_fine);

% unitypes = unique(celltype);
linewidth = 1;
% markers = {'+','o','*','x','s','d','^','v','>','<','p','h','.'};
% markers = {'o','*','s','d','v','>','<','p','h','.'};
% markers = {'o','s','s','o','o'};
markers = {'.','.','.','.'};
color = lines(length(ind_included_panel1));
color = color([2 1 3 4],:);

figure

g = subplot(2,2,1);
% plot(flipharm(pcs(1))*c_signal_pcastr.harmscr(:,pcs(1)),flipharm(pcs(2))*c_signal_pcastr.harmscr(:,pcs(2)),'.','Color',[.7 .7 .7]);
% set(gca,'XLim',[min(flipharm(pcs(1))*c_signal_pcastr.harmscr(:,pcs(1))) max(flipharm(pcs(1))*c_signal_pcastr.harmscr(:,pcs(1)))]*1.1,'YLim',[min(flipharm(pcs(2))*c_signal_pcastr.harmscr(:,pcs(2))) max(flipharm(pcs(2))*c_signal_pcastr.harmscr(:,pcs(2)))]*1.1)
hold on

xlabel(['PC ' num2str(pcs(1))])
ylabel(['PC ' num2str(pcs(2))])

legendstyles = nan(1,length(ind_included_panel1));

flipped_scores = repmat(flipharm,size(c_signal_pcastr.harmscr,1),1).*c_signal_pcastr.harmscr;

x_scores = (Rmat(pcs(1),:) * flipped_scores(:,1:3)')';
y_scores = (Rmat(pcs(2),:) * flipped_scores(:,1:3)')';
plot(x_scores,y_scores,'.','Color',[.7 .7 .7]);

fitcoef_rot = fitcoef;
fitcoef_rot(1:3,:) = Rmat*fitcoef_rot(1:3,:);

posFig = get(gcf,'Position');

for ilig = 1:length(ind_included_panel1)
    
    if sum(ind_included_panel1(ilig)==sites_for_harmonics)
        % Data used for harmonics
%                 legendstyles(ilig) = plot(c_signal_pcastr.harmscr(celltype(ind_harm) == unitypes(ilig),icol),c_signal_pcastr.harmscr(celltype(ind_harm) == unitypes(ilig),irow+1),'o','MarkerFaceColor',color(ilig,:),'LineWidth',linewidth);
        legendstyles(ilig) = plot(x_scores(celltype(ind_harm) == ind_included_panel1(ilig),:),y_scores(celltype(ind_harm) == ind_included_panel1(ilig),:),'.','Color',color(ilig,:));
    else
        % Data used for fitting
        legendstyles(ilig) = plot(fitcoef_rot(pcs(1),celltype(ind_fit) == ind_included_panel1(ilig)),fitcoef_rot(pcs(2),celltype(ind_fit) == ind_included_panel1(ilig)),markers{ilig},'MarkerEdgeColor',color(ilig,:),'LineWidth',linewidth);
    end
    
    set(gca,'XLim',[min(x_scores) max(x_scores)]*1.1)
    aspRatioFig = posFig(3)/posFig(4);
    posSubplot = get(gca,'Position');
    aspRatioSubplot = posSubplot(3)/posSubplot(4);
    set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)-.02)
        
end

input_names = {'No Lig', 'EGF', 'IGF', 'BTC'};
legend(g,legendstyles,input_names)

g = subplot(2,2,2);
% plot(flipharm(pcs(1))*c_signal_pcastr.harmscr(:,pcs(1)),flipharm(pcs(2))*c_signal_pcastr.harmscr(:,pcs(2)),'.','Color',[.7 .7 .7]);
% set(gca,'XLim',[min(flipharm(pcs(1))*c_signal_pcastr.harmscr(:,pcs(1))) max(flipharm(pcs(1))*c_signal_pcastr.harmscr(:,pcs(1)))]*1.1,'YLim',[min(flipharm(pcs(2))*c_signal_pcastr.harmscr(:,pcs(2))) max(flipharm(pcs(2))*c_signal_pcastr.harmscr(:,pcs(2)))]*1.1)
hold on

xlabel(['PC ' num2str(pcs(1))])
ylabel(['PC ' num2str(pcs(2))])

legendstyles = nan(1,length(ind_included_panel2));

plot(x_scores,y_scores,'.','Color',[.7 .7 .7]);

color = color([2 3 1 4],:);

for ilig = 1:length(ind_included_panel2)
    
    if sum(ind_included_panel2(ilig)==sites_for_harmonics)
        % Data used for harmonics
%                 legendstyles(ilig) = plot(c_signal_pcastr.harmscr(celltype(ind_harm) == unitypes(ilig),icol),c_signal_pcastr.harmscr(celltype(ind_harm) == unitypes(ilig),irow+1),'o','MarkerFaceColor',color(ilig,:),'LineWidth',linewidth);
        legendstyles(ilig) = plot(x_scores(celltype(ind_harm) == ind_included_panel2(ilig),:),y_scores(celltype(ind_harm) == ind_included_panel2(ilig),:),'.','Color',color(ilig,:));
    else
        % Data used for fitting
        legendstyles(ilig) = plot(fitcoef_rot(pcs(1),celltype(ind_fit) == ind_included_panel2(ilig)),fitcoef_rot(pcs(2),celltype(ind_fit) == ind_included_panel2(ilig)),markers{ilig},'MarkerEdgeColor',color(ilig,:),'LineWidth',linewidth);
    end
    
    set(gca,'XLim',[min(x_scores) max(x_scores)]*1.1)
    aspRatioFig = posFig(3)/posFig(4);
    posSubplot = get(gca,'Position');
    aspRatioSubplot = posSubplot(3)/posSubplot(4);
    set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot)-.02)
        
end

input_names = {'EGF alone', 'EGF + MEKi', 'EGF + AKTi'};
legend(g,legendstyles,input_names)

subplot(2,2,3)
hold on

color = lines(length(ind_included_panel1));
color = color([2 1 3 4],:);

for ilig = 1:length(ind_included_panel1)
    
    if sum(ind_included_panel1(ilig)==sites_for_harmonics)

        scores = (Rmat * flipped_scores(:,1:3)')';
        % Including PC1:
        pc_trajectories = scores(celltype(ind_harm) == ind_included_panel1(ilig),1:3)*(Rmat*harm_eval(:,1:3)');
%         pc_trajectories = scores(celltype(ind_harm) == ind_included_panel1(ilig),1:3)*(harm_eval(:,1:3)'); % not rotate harmonics!
        % Only PC2 & PC3:
%         pc_trajectories = scores(celltype(ind_harm) == ind_included_panel1(ilig),2:3)*(Rmat(2:3,:)*harm_eval(:,1:3)');
        mean_trajectories = mean(pc_trajectories,1);
        std_trajectories = std(pc_trajectories)./sqrt(size(pc_trajectories,1));
        
        tmpx = [times_fine'; flipud(times_fine')];
        tmpy = [mean_trajectories'+2*std_trajectories'; flipud(mean_trajectories'-2*std_trajectories')];
        
        plot(times_fine,mean_trajectories,'Color',color(ilig,:))
        
        ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
        set(ltmp, 'FaceColor', color(ilig,:)*0.1+0.9, 'EdgeColor', color(ilig,:)*0.1+0.9);
        ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
        set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', color(ilig,:)*0.3+0.7);
        
        plot(timestamp(range_ind),mean(c_signal(range_ind,celltype == ind_included_panel1(ilig)),2),'.','Color',color(ilig,:),'LineWidth',2)
        
    else
        % Data used for fitting
        legendstyles(ilig) = plot(fitcoef_rot(pcs(1),celltype(ind_fit) == ind_included_panel1(ilig)),fitcoef_rot(pcs(2),celltype(ind_fit) == ind_included_panel1(ilig)),markers{ilig},'MarkerEdgeColor',color(ilig,:),'LineWidth',linewidth);
    end
    
    xlabel('time')
%     set(gca,'XLim',time_range)

end

ylim = [-3 2]*1e-3;
% set(gca,'YLim',ylim)

subplot(2,2,4)
hold on

color = color([2 3 1 4],:);

for ilig = 1:length(ind_included_panel2)
    
    if sum(ind_included_panel2(ilig)==sites_for_harmonics)
    
        scores = (Rmat * flipped_scores(:,1:3)')';
        % Including PC1:
        pc_trajectories = scores(celltype(ind_harm) == ind_included_panel2(ilig),:)*(Rmat*harm_eval(:,1:3)');
%         pc_trajectories = scores(celltype(ind_harm) == ind_included_panel2(ilig),1:3)*(harm_eval(:,1:3)'); % not rotate harmonics!
        % Only PC2 & PC3:
%         pc_trajectories = scores(celltype(ind_harm) == ind_included_panel2(ilig),2:3)*(Rmat(2:3,:)*harm_eval(:,1:3)');
        mean_trajectories = mean(pc_trajectories,1);
        std_trajectories = std(pc_trajectories)./sqrt(size(pc_trajectories,1));
        
        tmpx = [times_fine'; flipud(times_fine')];
        tmpy = [mean_trajectories'+2*std_trajectories'; flipud(mean_trajectories'-2*std_trajectories')];
        
        plot(times_fine,mean_trajectories,'Color',color(ilig,:))
        
        ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
        set(ltmp, 'FaceColor', color(ilig,:)*0.1+0.9, 'EdgeColor', color(ilig,:)*0.1+0.9);
        ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
        set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', color(ilig,:)*0.3+0.7);
        
        plot(timestamp(range_ind),mean(c_signal(range_ind,celltype == ind_included_panel2(ilig)),2),'.','Color',color(ilig,:),'LineWidth',2)
        
    else
        % Data used for fitting
        scores = fitcoef_rot(1:3,:)';
        % Including PC1:
        pc_trajectories = scores(celltype(ind_fit) == ind_included_panel2(ilig),:)*(Rmat*harm_eval(:,1:3)');
%         pc_trajectories = scores(celltype(ind_fit) == ind_included_panel2(ilig),:)*(harm_eval(:,1:3)'); % not rotate harmonics!
        % Only PC2 & PC3:
%         pc_trajectories = scores(celltype(ind_fit) == ind_included_panel2(ilig),2:3)*(Rmat(2:3,:)*harm_eval(:,1:3)');
        mean_trajectories = mean(pc_trajectories,1);
        std_trajectories = std(pc_trajectories)./sqrt(size(pc_trajectories,1));
        
        tmpx = [times_fine'; flipud(times_fine')];
        tmpy = [mean_trajectories'+2*std_trajectories'; flipud(mean_trajectories'-2*std_trajectories')];
        
        plot(times_fine,mean_trajectories,'Color',color(ilig,:))
        
        ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
        set(ltmp, 'FaceColor', color(ilig,:)*0.1+0.9, 'EdgeColor', color(ilig,:)*0.1+0.9);
        ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
        set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', color(ilig,:)*0.3+0.7);
        
        plot(timestamp(range_ind),mean(c_signal(range_ind,celltype == ind_included_panel2(ilig)),2),'.','Color',color(ilig,:),'LineWidth',2)
    end
    
    xlabel('time')
%     set(gca,'XLim',time_range)

end

% set(gca,'YLim',ylim)