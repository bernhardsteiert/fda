%% Fit additional data-sets on harmonics generated from initial data-sets

close all
clear all
clc

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
addpath(grabdataPath)


% sites = [17 41 42 44 57];
% input_names = {'IGF','HGF-MEKi','HGF-AKTi','HGF','EPR'};
% sites_for_harmonics = [17 41 42 44 57];

sites = [1 2 4 10 17 57 64];
input_names = {'EGF-MEKi','EGF-AKTi','EGF','No Lig','IGF','EPR','BTC'};
sites_for_harmonics = [1 2 4 10 17 57 64];
% all ligands highest dose


times = cell(0);
signals = cell(0);
signals_raw = cell(0);
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
        signals_raw{end+1} = log10(intensity);
    else
        signals_raw{end+1} = intensity;
    end
    
    signals{end+1} = signals_raw{end} - repmat(nanmean(signals_raw{end},2),1,size(signals_raw{end},2));
    
    celltype = [celltype ones(1,size(intensity,2))*isite];
end

timestamp = times{1}; % same time sampling for all data sets
c_signal = cell2mat(signals);

return

% Exclude outliers
% HGF-AKTi, seventh data-set
exclude_site = [42];
exclude_signal = [7];
for iex = 1:length(exclude_site)
    tmpi = find(celltype == exclude_site(iex));
    ind_new = setdiff(1:size(c_signal,2),tmpi(1)+exclude_signal(iex)-1);
    c_signal = c_signal(:,ind_new);
    celltype = celltype(ind_new);
end

return

%% Plot raw data of defined ligands
close all
figure

c_signal_raw = cell2mat(signals_raw);

plot_ligs = sites;
plot_name = input_names;
nrows = 3;
ncols = 3;

nplot = 20;

for iplot = 1:length(plot_ligs)
    subplot(nrows,ncols,iplot)
    
    plot(repmat(timestamp,1,sum(celltype == plot_ligs(iplot))),c_signal_raw(:,celltype == plot_ligs(iplot)),'g','color',[0.7 0.7 0.7])
    colored_ind = find(celltype == plot_ligs(iplot));
    hold on
%     plot(repmat(timestamp,1,nplot),c_signal_raw(:,colored_ind(1:nplot)))
    plot(timestamp,nanmean(c_signal_raw(:,celltype == plot_ligs(iplot)),2),'color','k','LineWidth',2)
    title(plot_name{iplot})
    
    ylim = [-1 1]*.04;
    if ~log_trafo
        ylim = 10.^ylim;
    end
    set(gca,'XLim',[0 650],'YLim',ylim)
    set(gca,'XLim',[0 650])

end

%% Plot every data set with distinct color
close all

% plot_sites = site;
plot_sites = sites;

rowstocols = 1;
nrows = ceil(length(plot_sites)^rowstocols);
ncols = ceil(length(plot_sites) / nrows);

figure

for ip = 1:length(plot_sites)
    subplot(nrows,ncols,ip)
    
    c_signal_single = c_signal(:,celltype == plot_sites(ip));
    
    first_n = 100; % Plot only first_n data-sets
    
    first_n = min(first_n,size(c_signal_single,2));
%     f = figure;
%     set(f,'DefaultAxesColorOrder',jet(first_n))
    
    plot(timestamp,c_signal_single(:,1:first_n))
%     title(['Site ' num2str(plot_sites(ip))])
    title(input_names{ip})
    
    set(gca,'XLim',[200 650])
    
    hold on
    plot(timestamp,nanmean(c_signal_single,2),'--k')
    
    plot([120 120],[-0.06 0.06],'b--')
    set(gca,'YLim',[-0.06 0.06])
    
%     if length(plot_sites) > 1
%         waitforbuttonpress;
%         close gcf
%     end
end

%% Generate spline fits to data-sets given in sites_for_harmonics
close all

nbasis = 40;
% time_range = [min(timestamp) max(timestamp)];
time_range = [200 650];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

ind_harm = ismember(celltype,sites_for_harmonics);
ind_fit = ~ind_harm;

basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
smoothed_data = smooth_basis(timestamp(range_ind),c_signal(range_ind,ind_harm),basis);

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal(1,ind_harm),2)))
hold on

plot(smoothed_data)
plot(timestamp(range_ind),c_signal(range_ind,ind_harm),'o')

%% Make FPCA with data generated in previous block
close all

nharm = 3;
c_signal_pcastr = pca_fd(smoothed_data, nharm);
% c_signal_pcastr = pca_fd(smoothed_data, nharm, fdPar(basis, int2Lfd(2), 0), 0); % WITHOUT CENTERING!!

plot_pca_fd(c_signal_pcastr, 1, 0)
% plot(c_signal_pcastr.meanfd)

% c_signal_rotpcastr = varmx_pca(c_signal_pcastr);
% plot_pca_fd(c_signal_rotpcastr, 1, 0)

%% Plot: Eigenfunctions
close all

rowstocols = 0.5;
nrows = ceil(nharm^rowstocols);
ncols = ceil(nharm / nrows);

time_range = [200 650];

flipharm = ones(1,nharm);
% flipharm(1:8) = [1 -1 1 -1 -1 1 1 1];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;
times_fine = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),501);

harm_eval = 2 * repmat(sqrt(c_signal_pcastr.values(1:nharm))'.*flipharm(1:nharm),length(times_fine),1) .* eval_fd(c_signal_pcastr.harmfd,times_fine);

for iplot = 1:nharm
    subplot(nrows,ncols,iplot)
    
    plot(times_fine,harm_eval(:,iplot))
    xlabel(['Harmonic ' num2str(iplot)])
    set(gca,'XLim',time_range)
    set(gca,'YLim',[min(min(harm_eval)) max(max(harm_eval))])
    
    hold on
    plot(time_range,[0 0],'--')
end

%% Remove all nharm from signal --> only stochastic oscillations remain
close all

c_signal_woNharm = c_signal(range_ind,ind_harm)-eval_fd(c_signal_pcastr.fdhatfd,timestamp(range_ind));
celltypeharm = celltype(ind_harm);

plot_sites = sites;

rowstocols = 1;
nrows = ceil(length(plot_sites)^rowstocols);
ncols = ceil(length(plot_sites) / nrows);

figure

for ip = 1:length(plot_sites)
    subplot(nrows,ncols,ip)
    
    c_signal_single = c_signal_woNharm(:,celltypeharm == plot_sites(ip));
    
    first_n = 100;
    
    first_n = min(first_n,size(c_signal_single,2));
    
    plot(timestamp(range_ind),c_signal_single(:,1:first_n))
    title(input_names{ip})
    
    set(gca,'XLim',[200 650])
    
    hold on
    plot(timestamp(range_ind),nanmean(c_signal_single,2),'--k')
    
    plot([120 120],[-0.04 0.04],'b--')
%     set(gca,'YLim',[-0.04 0.04])
end

%% Generate spline fits to data-sets given in sites_for_harmonics (for remaining variation)
close all

nbasis = 40;

% ind_harm = ismember(celltype,sites_for_harmonics);
% ind_fit = ~ind_harm;

smoothed_data_woNharm = smooth_basis(timestamp(range_ind),c_signal_woNharm,basis);

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal(1,ind_harm),2)))
hold on

plot(smoothed_data_woNharm)
plot(timestamp(range_ind),c_signal_woNharm,'o')

%% Plot: Histogramm of distance to origin
close all

rad_dist_thres = 0.03;

figure
hold on

rowstocols = 0.6;
nrows = ceil(length(plot_sites)^rowstocols);
ncols = ceil(length(plot_sites) / nrows);

radial_dist = sqrt(sum(getcoef(smoothed_data_woNharm).^2,1));

posFig = get(gcf,'Position');
% posFig(3) = posFig(3)/2;
% posFig(4) = posFig(4)*2;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./15);

for ip = 1:length(plot_sites)
    subplot(nrows,ncols,ip)
    
    baredges = linspace(0,max(radial_dist)+.01,21);
    bar(baredges,histc(radial_dist(celltypeharm == plot_sites(ip)),baredges));
    
    title(input_names{ip})
    
    set(gca,'XLim',[0 max(radial_dist)+.01])
    if ip == length(plot_sites)
        xlabel('radial distance')
    end
    if ip == 1
        ylabel('absolute frequency')
    end
    
    hold on
    
%     plot([rad_dist_thres rad_dist_thres],get(gca,'YLim'),'--')
end

%% Plot: Histogramm of distance to origin (new - overlayed)
close all

groups = {[1 2 4 10], [4 10 17 57 64]};
resort = {[2 3 1 4], [1 4 2 3 5]};

figure
hold on

rowstocols = 0.6;
nrows = ceil(length(groups)^rowstocols);
ncols = ceil(length(groups) / nrows);

radial_dist = sqrt(sum(getcoef(smoothed_data_woNharm).^2,1));

for ig = 1:length(groups)
    g = subplot(nrows,ncols,ig);
    hold on
    
    mygroup = groups{ig};
    
    legend_names = {};
    legendstyles = nan(1,length(mygroup));
    color = lines(length(mygroup));
    color(4,:) = 0;
    for ip = 1:length(mygroup)

%         baredges = linspace(0,max(radial_dist)+.01,26);
%         barheight = histc(radial_dist(celltypeharm == mygroup(ip)),baredges)./sum(celltypeharm == mygroup(ip));
%         legendstyles(ip) = plot(baredges,barheight,'Color',color(resort{ig}(ip),:));

        [f,xi] = ksdensity(radial_dist(celltypeharm == mygroup(ip)));
        legendstyles(ip) = plot(xi,f,'Color',color(resort{ig}(ip),:));

        legend_names{end+1} = (input_names{mygroup(ip) == sites_for_harmonics});

    end
    
    legend(g,legendstyles,legend_names)
    
    set(gca,'XLim',[0 max(radial_dist)+.01])
    set(gca,'YLim',[0 200])

    xlabel('radial distance')
%     ylabel('relative frequency')
    ylabel('estimated density')

end

%% Plot traces with gray level corrsesponding to radial distance
close all

rowstocols = 0.6;
nrows = ceil(length(plot_sites)^rowstocols);
ncols = ceil(length(plot_sites) / nrows);

figure

posFig = get(gcf,'Position');
% posFig(3) = posFig(3)/2;
% posFig(4) = posFig(4)*2;
set(gcf,'Position',posFig)
set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./15);

ncolor = 201;
color = repmat(linspace(0,1,ncolor),3,1);

color = color(:,end:-1:1); % Gray scale - darkness depending on score

% TODO: As Pat suggested: Plot 4 representing lines / subplot in red - saturation depending on score

% color = (abs(1-color .* lines(ncolor)')).^2; % Colored - saturation depending on score
radial_space = linspace(min(radial_dist),max(radial_dist),ncolor);
[radial_dist_sorted ind_sort_radial] = sort(radial_dist);

for ip = 1:length(plot_sites)
    subplot(nrows,ncols,ip)
    hold on
    
    c_signal_single = c_signal_raw(:,ind_sort_radial);
    c_signal_single = c_signal_single(:,(celltypeharm(ind_sort_radial) == plot_sites(ip)));
    radial_dist_single = radial_dist_sorted(celltypeharm(ind_sort_radial) == plot_sites(ip));
    
    for i = 1:size(c_signal_single,2)
        [tmp color_ind] = min(abs(radial_dist_single(i) - radial_space));
        plot(timestamp,c_signal_single(:,i),'Color',color(:,color_ind))
    end
    
    plot(get(gca,'XLim'),[0 0],'--k')
    
    title(input_names{ip})
    
    set(gca,'XLim',[50 650])
    
    plot([120 120],[-0.04 0.04],'b--')
%     set(gca,'YLim',[-0.04 0.04])
end


%% Plot traces under/over a defined radial distance threshold seperately
close all

rowstocols = 1;
nrows = ceil(length(plot_sites)^rowstocols);
ncols = ceil(length(plot_sites) / nrows)*2;

figure

for ip = 1:length(plot_sites)
    subplot(nrows,ncols,(ip-1)*2+1)
    
%     c_signal_single = c_signal_woNharm(:,(celltypeharm == plot_sites(ip)) & (radial_dist <= rad_dist_thres));
    c_signal_single = c_signal_raw(:,(celltypeharm == plot_sites(ip)) & (radial_dist <= rad_dist_thres));

    first_n = 10;
    
    first_n = min(first_n,size(c_signal_single,2));
    
    try
        plot(timestamp,c_signal_single(:,1:first_n))
    end
    
    hold on
    plot(get(gca,'XLim'),[0 0],'--k')
    
    title(input_names{ip})
    
    set(gca,'XLim',[50 650])
    
    plot([120 120],[-0.04 0.04],'b--')
%     set(gca,'YLim',[-0.04 0.04])
end

for ip = 1:length(plot_sites)
    subplot(nrows,ncols,ip*2)
    
%     c_signal_single = c_signal_woNharm(:,(celltypeharm == plot_sites(ip)) & (radial_dist > rad_dist_thres));
    c_signal_single = c_signal_raw(:,(celltypeharm == plot_sites(ip)) & (radial_dist > rad_dist_thres));

    first_n = 10;
    
    first_n = min(first_n,size(c_signal_single,2));
    
    try
        plot(timestamp,c_signal_single(:,1:first_n))
    end
    
    hold on
    plot(get(gca,'XLim'),[0 0],'--k')
    
    title(input_names{ip})
    
    set(gca,'XLim',[50 650])
    
    plot([120 120],[-0.04 0.04],'b--')
%     set(gca,'YLim',[-0.04 0.04])
end


%% Plot: %variance explained vs. #basis functions
close all

thres_var = 0.9;

cumprobs = cumsum([0;c_signal_pcastr.varprop]);
[tmp thres_ind] = min(abs(cumprobs - thres_var));
if cumprobs(thres_ind)-thres_var < 0
    thres_ind = thres_ind + 1;
end

plot(0:length(c_signal_pcastr.varprop),cumprobs)
hold on
plot(0:thres_ind-1,ones(1,thres_ind)*thres_var,'--')
plot([thres_ind-1 thres_ind-1],[0 cumprobs(thres_ind)],'--')

xlabel('fPCA basis functions')
ylabel('cumulative variance explained')

fprintf('To explain at least %s variance, use %i fPCA basis functions.\n\n',num2str(thres_var,3),thres_ind-1);

%% Plot: Eigenvalues linear trend
close all

linearrange = 8:15;

axb = polyfit(linearrange,log10(c_signal_pcastr.values(linearrange))',1);
plot(1:max(linearrange),log10(c_signal_pcastr.values(1:max(linearrange))),'ko-')
hold on
plot(1:max(linearrange),(1:max(linearrange))*axb(1) + axb(2),'k--')

xlabel('Eigenvalue Number')
ylabel('log10(Eigenvalue)')

%% Fit additional data with basis from fPCA
close all

harm_basis = create_fd_basis(c_signal_pcastr.harmfd);
mean_fit = eval_fd(c_signal_pcastr.meanfd,timestamp(range_ind));
smoothed_additional = smooth_basis(timestamp(range_ind),c_signal(range_ind,ind_fit)-repmat(mean_fit,1,sum(ind_fit)),harm_basis);

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal(1,ind_fit),2)))
hold on

plot(smoothed_additional+c_signal_pcastr.meanfd)
plot(timestamp(range_ind),c_signal(range_ind,ind_fit),'o')

%% Check if data is fitted by harmonics
close all

nperplot = 8;
nsubplots = ceil(size(c_signal(1,ind_fit),2)./nperplot);
rowstocols = 0.5;
time_range = [50 650];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;
times_fine = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),501);
harm_fine = eval_basis(harm_basis,times_fine);

nrows = ceil(nsubplots^rowstocols);
ncols = ceil(nsubplots / nrows);

fitcoef = getcoef(smoothed_additional);
data_fpca_repr = fitcoef'*harm_fine';
mean_fine = eval_fd(c_signal_pcastr.meanfd,times_fine);

ind_fit_no = find(ind_fit);

for iplot = 1:nsubplots
    subplot(nrows,ncols,iplot)
    
    inds = ((iplot-1)*nperplot+1):min([size(c_signal(1,ind_fit),2) iplot*nperplot]);
    plot(times_fine,repmat(mean_fine,1,length(inds))+data_fpca_repr(inds,:)')
    hold on
    plot(timestamp(range_ind),c_signal(range_ind,ind_fit_no(inds)),'o')

end

%% Plot: Triagonal Matrix of PCs (the bold points are fitted)
close all

max_pc = 3;

% If no harmonics were fitted:
fitcoef = zeros(size(c_signal_pcastr.harmscr'));

figure
color = lines(length(signals));
legendstyles = nan(1,length(signals));
hold on

unitypes = unique(celltype);
linewidth = 1;
% markers = {'+','o','*','x','s','d','^','v','>','<','p','h','.'};
markers = {'o','*','s','d','v','>','<','p','h','.'};
% markers = {'o','s','s','o','o'};

xpos = .07;
ypos = .04;

for irow = 1:max_pc-1
    for icol = 1:irow
        
        h = subplot(max_pc-1,max_pc-1,(irow-1)*(max_pc-1)+icol);
        pos = get(h,'Pos');
        set(h,'Pos',[pos(1)-xpos*(max_pc-1-icol)/max_pc pos(2)-ypos*irow/max_pc pos(3)*1.2 pos(4)*1.2])
        hold on
        if irow == max_pc-1;
            xlabel(['PC ' num2str(icol)])
        else
            set(gca,'XTickLabel',[])
        end
        if icol == 1
            ylabel(['PC ' num2str(irow+1)])
        else
            set(gca,'YTickLabel',[])
        end
        
        for ilig = 1:length(signals)
            if sum(unitypes(ilig)==sites_for_harmonics)
                % Data used for harmonics
%                 legendstyles(ilig) = plot(c_signal_pcastr.harmscr(celltype(ind_harm) == unitypes(ilig),icol),c_signal_pcastr.harmscr(celltype(ind_harm) == unitypes(ilig),irow+1),'o','MarkerFaceColor',color(ilig,:),'LineWidth',linewidth);
                legendstyles(ilig) = plot(c_signal_pcastr.harmscr(celltype(ind_harm) == unitypes(ilig),icol),c_signal_pcastr.harmscr(celltype(ind_harm) == unitypes(ilig),irow+1),markers{ilig},'MarkerEdgeColor',color(ilig,:),'LineWidth',linewidth+1);
            else
                % Data used for fitting
                legendstyles(ilig) = plot(fitcoef(icol,celltype(ind_fit) == unitypes(ilig)),fitcoef(irow+1,celltype(ind_fit) == unitypes(ilig)),'s','MarkerEdgeColor',color(ilig,:),'LineWidth',linewidth+1);
            end
            
        end
        
        xrange = max([c_signal_pcastr.harmscr(:,icol);fitcoef(icol,:)']) - min([c_signal_pcastr.harmscr(:,icol);fitcoef(icol,:)']);
        yrange = max([c_signal_pcastr.harmscr(:,irow+1);fitcoef(irow+1,:)']) - min([c_signal_pcastr.harmscr(:,irow+1);fitcoef(irow+1,:)']);
        scalefac = .05;
        
        set(gca,'XLim',[min([c_signal_pcastr.harmscr(:,icol);fitcoef(icol,:)']) - xrange*scalefac max([c_signal_pcastr.harmscr(:,icol);fitcoef(icol,:)']) + xrange*scalefac])
        set(gca,'YLim',[min([c_signal_pcastr.harmscr(:,irow+1);fitcoef(irow+1,:)']) - yrange*scalefac max([c_signal_pcastr.harmscr(:,irow+1);fitcoef(irow+1,:)']) + yrange*scalefac])
        
    end
end

g = subplot(max_pc-1,max_pc-1,max_pc-1);
set(gca,'Visible','off')
legend(g,legendstyles,input_names)
