%% Fit additional data-sets on harmonics generated from initial data-sets

close all
clear all
clc

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
addpath(grabdataPath)

% Get properties of sites by calling siteprop(site)
% sites = [1 2 4:9 17:-1:12 57:-1:52];
% sites_for_harmonics = [4:9 17:-1:12 57:-1:52];

% All ligands
% sites = [4:9 17:-1:12 24:29 37:-1:32 44:49 57:-1:52 64:69];
% sites_for_harmonics = [4:9 17:-1:12 24:29 37:-1:32 44:49 57:-1:52 64:69];

% All ligands + EGF with Inhibitors
sites = [1 2 4:9 17:-1:12 24:29 37:-1:32 44:49 57:-1:52 64:69];
sites_for_harmonics = [4:9 17:-1:12 24:29 37:-1:32 44:49 57:-1:52 64:69];


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

%% Generate spline fits to data-sets given in sites_for_harmonics
close all

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

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal(1,ind_harm),2)))
hold on

plot(smoothed_data)
plot(timestamp(range_ind),c_signal(range_ind,ind_harm),'o')

%% Make FPCA with data generated in previous block - wip
close all

nharm = 8;
% c_signal_pcastr = pca_fd(smoothed_data, nharm);
c_signal_pcastr = pca_fd(smoothed_data, nharm, fdPar(basis, int2Lfd(2), 0), 0); % WITHOUT CENTERING!!

plot_pca_fd(c_signal_pcastr, 1, 0)

% c_signal_rotpcastr = varmx_pca(c_signal_pcastr);
% plot_pca_fd(c_signal_rotpcastr, 1, 0)

%% Plot: Eigenfunctions
close all

rowstocols = 0.5;
nrows = ceil(nharm^rowstocols);
ncols = ceil(nharm / nrows);

time_range = [50 200];

flipharm = ones(1,nharm);
% flipharm(1:4) = [1 1 -1 -1];

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


%% Plot: Reproduce figure from poster with colored ligand level
close all

% Define principal components to be plotted
pcs = [1 2];

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
color_doses = 10.^linspace(log10(lig_min),log10(lig_max),ncolor);

rowstocols = 0.5;
nrows = ceil((length(uni_lig)+1)^rowstocols);
ncols = ceil((length(uni_lig)+1) / nrows);

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
    
    plot(c_signal_pcastr.harmscr(:,pcs(1)),c_signal_pcastr.harmscr(:,pcs(2)),'.','Color',[.7 .7 .7]);
    
    for isite = tmpind
        % Colored
        [tmp color_ind] = min(abs(color_doses-site_lig_dose(sites_remain(isite))));
        mycolor = colmap(color_ind,:);
        plot(c_signal_pcastr.harmscr(celltype(ind_harm) == sites_for_harmonics(sites_remain(isite)),pcs(1)),c_signal_pcastr.harmscr(celltype(ind_harm) == sites_for_harmonics(sites_remain(isite)),pcs(2)),'.','Color',mycolor);
    end
    
    set(gca,'XLim',[min(c_signal_pcastr.harmscr(:,pcs(1))) max(c_signal_pcastr.harmscr(:,pcs(1)))]*1.1,'YLim',[min(c_signal_pcastr.harmscr(:,pcs(2))) max(c_signal_pcastr.harmscr(:,pcs(2)))]*1.1)
    
end

h = get(gca);
subplot(nrows,ncols,ilig+1)
set(gca,'CLim',log10([lig_min lig_max]))
colormap(flipud(jet(ncolor)))
colorbar('Location','North','XTick',log10([2.5 5 10 20 50 100]),'XTickLabel',[2.5 5 10 20 50 100])
set(gca,'Visible','off')


%% Plot: Show effect of different eigenfunctions by color circle
close all

circ_thres = .9;

subplot(4,4,[1,2,5,6])
plot(c_signal_pcastr.harmscr(:,pcs(1)),c_signal_pcastr.harmscr(:,pcs(2)),'.','Color',[.7 .7 .7]);

set(gca,'XLim',[min(c_signal_pcastr.harmscr(:,pcs(1))) max(c_signal_pcastr.harmscr(:,pcs(1)))]*1.1,'YLim',[min(c_signal_pcastr.harmscr(:,pcs(2))) max(c_signal_pcastr.harmscr(:,pcs(2)))]*1.1)

ylabel(['PC ' num2str(pcs(2))])

hold on

% [dists dist_ind] = sort((c_signal_pcastr.harmscr(:,pcs(1))-mean(c_signal_pcastr.harmscr(:,pcs(1)))).^2 + (c_signal_pcastr.harmscr(:,pcs(2))-mean(c_signal_pcastr.harmscr(:,pcs(2)))).^2);
[dists dist_ind] = sort((1./c_signal_pcastr.values(pcs(1))*(c_signal_pcastr.harmscr(:,pcs(1))-mean(c_signal_pcastr.harmscr(:,pcs(1))))).^2 + (1./c_signal_pcastr.values(pcs(2))*(c_signal_pcastr.harmscr(:,pcs(2))-mean(c_signal_pcastr.harmscr(:,pcs(2))))).^2);

tcirc = 0:0.01:1;
plot(mean(c_signal_pcastr.harmscr(:,pcs(1)))+c_signal_pcastr.values(pcs(1))*sqrt(dists(floor(length(dists)*circ_thres)))*sin(2*pi*tcirc),mean(c_signal_pcastr.harmscr(:,pcs(2)))+c_signal_pcastr.values(pcs(2))*sqrt(dists(floor(length(dists)*circ_thres)))*cos(2*pi*tcirc),'k--')


subplot(4,4,[3,4,7,8])
plot(c_signal_pcastr.harmscr(dist_ind(floor(length(dists)*circ_thres):end),pcs(1)),c_signal_pcastr.harmscr(dist_ind(floor(length(dists)*circ_thres):end),pcs(2)),'.','Color',[.7 .7 .7]);

set(gca,'XLim',[min(c_signal_pcastr.harmscr(:,pcs(1))) max(c_signal_pcastr.harmscr(:,pcs(1)))]*1.1,'YLim',[min(c_signal_pcastr.harmscr(:,pcs(2))) max(c_signal_pcastr.harmscr(:,pcs(2)))]*1.1)

xlabel(['PC ' num2str(pcs(1))])

hold on

radius = linspace(0,sqrt(dists(floor(length(dists)*circ_thres))),6);
angle = 0:0.05:1;

X = c_signal_pcastr.values(pcs(1))*radius'*cos(2*pi*angle);
Y = c_signal_pcastr.values(pcs(2))*radius'*sin(2*pi*angle);

H = repmat(linspace(0,1,length(angle)),length(radius),1);
S = repmat(linspace(0,1,length(radius))',1,length(angle));
V = .9*ones(length(radius),length(angle));

hsvImage = cat(3,H,S,V);
C = hsv2rgb(hsvImage);

for r = 1:length(radius)-1
    for a = 1:length(angle)-1
        x = X(r:r+1,a:a+1);
        y = Y(r:r+1,a:a+1);
        patch(mean(c_signal_pcastr.harmscr(:,pcs(1)))+x([1 2 4 3]),mean(c_signal_pcastr.harmscr(:,pcs(2)))+y([1 2 4 3]),C(r,a,:))
    end
end

subplot(4,4,[9,10,11,13,14,15])

hold on

for r = 2:length(radius)-1
    for a = 2:length(angle)-1
        x = X(r:r+1,a:a+1);
        y = Y(r:r+1,a:a+1);
        plot(times_fine,(mean(c_signal_pcastr.harmscr(:,pcs(1)))+mean(x(:)))*harm_eval(:,pcs(1))+(mean(c_signal_pcastr.harmscr(:,pcs(2)))+mean(y(:)))*harm_eval(:,pcs(2)),'Color',C(r,a,:))
    end
end

plot(times_fine,mean(c_signal_pcastr.harmscr(:,pcs(1)))*harm_eval(:,pcs(1))+mean(c_signal_pcastr.harmscr(:,pcs(2)))*harm_eval(:,pcs(2)),'k','LineWidth',2)

xlabel('time')
ylabel('reconstructed signal')

subplot(4,4,12)

plot(times_fine,harm_eval(:,pcs(1)))
title(['Harmonic ' num2str(pcs(1))])
set(gca,'XLim',time_range)
set(gca,'YLim',[min(min(harm_eval)) max(max(harm_eval))]*1.1)
set(gca,'XTick',[],'YTick',[])

hold on
plot(time_range,[0 0],'--')

subplot(4,4,16)

plot(times_fine,harm_eval(:,pcs(2)))
title(['Harmonic ' num2str(pcs(2))])
set(gca,'XLim',time_range)
set(gca,'YLim',[min(min(harm_eval)) max(max(harm_eval))]*1.1)
set(gca,'XTick',[],'YTick',[])

hold on
plot(time_range,[0 0],'--')

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
time_range = [50 200];

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
% figure

max_pc = 4;

% If no harmonics were fitted:
% fitcoef = zeros(size(c_signal_pcastr.harmscr'));

figure
color = lines(length(signals));
% legendstyles = nan(1,length(signals));
legendstyles = nan(1,3);
hold on

unitypes = unique(celltype);
linewidth = 1;
% markers = {'+','o','*','x','s','d','^','v','>','<','p','h','.'};
% markers = {'o','*','s','d','v','>','<','p','h','.'};
% markers = {'o','s','s','o','o'};
markers = {'o','o','s','s','s','s','s','s','s','s','s','s','s','s','s','s','s','s','s','s'};
color = lines(3);

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
        
        cinh = 1;
        
        resort = [];
        for ilig = 1:length(signals)
            if unitypes(ilig) == 4 % EGF Alone High dose
                resort = [resort ilig];
            elseif sum(unitypes(ilig)==sites_for_harmonics)
                % Data used for harmonics
                resort = [ilig resort];
            else
                % Data used for fitting
                resort = [resort ilig];
            end
            
        end
        
        for ilig = resort
            if unitypes(ilig) == 4 % EGF Alone High dose
                legendstyles(cinh) = plot(c_signal_pcastr.harmscr(celltype(ind_harm) == unitypes(ilig),icol),c_signal_pcastr.harmscr(celltype(ind_harm) == unitypes(ilig),irow+1),'.','Color',color(cinh,:),'LineWidth',linewidth);
            elseif sum(unitypes(ilig)==sites_for_harmonics)
                % Data used for harmonics
%                 legendstyles(ilig) = plot(c_signal_pcastr.harmscr(celltype(ind_harm) == unitypes(ilig),icol),c_signal_pcastr.harmscr(celltype(ind_harm) == unitypes(ilig),irow+1),'o','MarkerFaceColor',color(ilig,:),'LineWidth',linewidth);
                plot(c_signal_pcastr.harmscr(celltype(ind_harm) == unitypes(ilig),icol),c_signal_pcastr.harmscr(celltype(ind_harm) == unitypes(ilig),irow+1),'.','Color',[.7 .7 .7],'LineWidth',linewidth);
            else
                % Data used for fitting
                legendstyles(cinh) = plot(fitcoef(icol,celltype(ind_fit) == unitypes(ilig)),fitcoef(irow+1,celltype(ind_fit) == unitypes(ilig)),'.','MarkerEdgeColor',color(cinh,:),'LineWidth',linewidth);
                cinh = cinh + 1;
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
% legend(g,legendstyles,input_names)

leg_str = cell(0);
for isite = [1 2 4]
    s = siteprop(isite);
    leg_str{end+1} = [s.lig_name num2str(s.lig_dose) s.inh_name];
end

legend(g,legendstyles,leg_str)


%% Mean-subtract data for capturing stochastic fluctuations
meansub_signals = signals;

for isite = 1:length(sites)
    meansub_signals{isite} = signals{isite} - repmat(nanmean(signals{isite},2),1,size(signals{isite},2));
end

c_signal_meansub = cell2mat(signals);


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
smoothed_data = smooth_basis(timestamp(range_ind),c_signal_meansub(range_ind,ind_harm),basis);

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal_meansub(1,ind_harm),2)))
hold on

plot(smoothed_data)
plot(timestamp(range_ind),c_signal_meansub(range_ind,ind_harm),'o')

%% Make FPCA with data generated in previous block
close all

nharm = 3;
c_signal_meansub_pcastr = pca_fd(smoothed_data, nharm);

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

harm_eval = 2 * repmat(sqrt(c_signal_meansub_pcastr.values(1:nharm))'.*flipharm(1:nharm),length(times_fine),1) .* eval_fd(c_signal_meansub_pcastr.harmfd,times_fine);

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

c_signal_meansub_woNharm = c_signal_meansub(range_ind,ind_harm)-eval_fd(c_signal_meansub_pcastr.fdhatfd,timestamp(range_ind));
celltypeharm = celltype(ind_harm);

ncols = 6;
nrows = ceil(length(sites_remain)/ncols);

figure

for ilig = 1:length(sites_remain)
    ip = sites_for_harmonics(ilig);
    subplot(nrows,ncols,ilig)
    
    c_signal_meansub_single = c_signal_meansub_woNharm(:,celltypeharm == ip);
    
    first_n = 100;
    
    first_n = min(first_n,size(c_signal_meansub_single,2));
    
    plot(timestamp(range_ind),c_signal_meansub_single(:,1:first_n))
    title([site_lig_name{ilig} num2str(site_lig_dose(ilig))])
    
    set(gca,'XLim',[200 650])
    
    hold on
    plot(timestamp(range_ind),nanmean(c_signal_meansub_single,2),'--k')
    
    plot([120 120],[-0.04 0.04],'b--')
    set(gca,'YLim',[-0.04 0.04])
end

%% Generate spline fits to data-sets given in sites_for_harmonics (for remaining variation)
close all

nbasis = 40;

% ind_harm = ismember(celltype,sites_for_harmonics);
% ind_fit = ~ind_harm;

smoothed_data_woNharm = smooth_basis(timestamp(range_ind),c_signal_meansub_woNharm,basis);

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal_meansub(1,ind_harm),2)))
hold on

plot(smoothed_data_woNharm)
plot(timestamp(range_ind),c_signal_meansub_woNharm,'o')

%% Plot: Histogramm of distance to origin
close all

rad_dist_thres = 0.025;

radial_dist = sqrt(sum(getcoef(smoothed_data_woNharm).^2,1));

ncols = 6;
nrows = ceil(length(sites_remain)/ncols);

figure
hold on

for ilig = 1:length(sites_remain)
    subplot(nrows,ncols,ilig)
    
    baredges = linspace(0,max(radial_dist)+.01,21);
    bar(baredges,histc(radial_dist(celltypeharm == sites_for_harmonics(ilig)),baredges));
    
    title([site_lig_name{ilig} num2str(site_lig_dose(ilig))])
    
    set(gca,'XLim',[0 max(radial_dist)+.01])
    if ip == length(sites_remain)
        xlabel('radial distance')
    end
    if ip == ceil(length(sites_remain)/2)
        ylabel('absolute frequency')
    end
    
    hold on
    
    plot([rad_dist_thres rad_dist_thres],get(gca,'YLim'),'--')
end

%% Plot traces under/over a defined radial distance threshold seperately
close all

ncols = 6;
nrows = ceil(length(sites_remain)/ncols);

figure
hold on

% for ip = 1:length(plot_sites)
%     subplot(nrows,ncols,(ip-1)*2+1)
%     
%     c_signal_meansub_single = c_signal_meansub_woNharm(:,(celltypeharm == plot_sites(ip)) & (radial_dist <= rad_dist_thres));
%     
%     first_n = 20;
%     
%     first_n = min(first_n,size(c_signal_meansub_single,2));
%     
%     try
%         plot(timestamp(range_ind),c_signal_meansub_single(:,1:first_n))
%     end
%     
%     hold on
%     plot(get(gca,'XLim'),[0 0],'--k')
%     
%     title(input_names{ip})
%     
%     set(gca,'XLim',[200 650])
%     
%     plot([120 120],[-0.04 0.04],'b--')
%     set(gca,'YLim',[-0.04 0.04])
% end

for ilig = 1:length(sites_remain)
    subplot(nrows,ncols,ilig)
    
    c_signal_meansub_single = c_signal_meansub_woNharm(:,(celltypeharm == sites_for_harmonics(ilig)) & (radial_dist > rad_dist_thres));
    
    first_n = 20;
    
    first_n = min(first_n,size(c_signal_meansub_single,2));
    
    try
        plot(timestamp(range_ind),c_signal_meansub_single(:,1:first_n))
    end
    
    hold on
    plot(get(gca,'XLim'),[0 0],'--k')
    
    title([site_lig_name{ilig} num2str(site_lig_dose(ilig))])
    
    set(gca,'XLim',[200 650])
    
    plot([120 120],[-0.04 0.04],'b--')
    set(gca,'YLim',[-0.04 0.04])
end

%% Make comprehensive PCA with scores of fPCA and radial distance combined
% Comment: This might not be a good idea ... maybe better as Pat already
% suggested: Take radial axis as another dimension!
close all

% When doing PCA, scores are compared to radial distance ... as these are
% not directly comparable, scale influence of oscillation accordingly
rad_scale_fac = 8;

pcainput = [c_signal_pcastr.harmscr(:,1:4) rad_scale_fac*radial_dist'];

[pcacoeff pcascores] = princomp(pcainput);

% Define principal components to be plotted
pcs = [1 2];

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
color_doses = 10.^linspace(log10(lig_min),log10(lig_max),ncolor);

rowstocols = 0.5;
nrows = ceil((length(uni_lig)+1)^rowstocols);
ncols = ceil((length(uni_lig)+1) / nrows);

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
    
    plot(pcascores(:,pcs(1)),pcascores(:,pcs(2)),'.','Color',[.7 .7 .7]);
    
    for isite = tmpind
        % Colored
        [tmp color_ind] = min(abs(color_doses-site_lig_dose(sites_remain(isite))));
        mycolor = colmap(color_ind,:);
        plot(pcascores(celltype(ind_harm) == sites_for_harmonics(sites_remain(isite)),pcs(1)),pcascores(celltype(ind_harm) == sites_for_harmonics(sites_remain(isite)),pcs(2)),'.','Color',mycolor);
    end
    
    set(gca,'XLim',[min(pcascores(:,pcs(1))) max(pcascores(:,pcs(1)))]*1.1,'YLim',[min(pcascores(:,pcs(2))) max(pcascores(:,pcs(2)))]*1.1)
    
end

h = get(gca);
subplot(nrows,ncols,ilig+1)
set(gca,'CLim',log10([lig_min lig_max]))
colormap(flipud(jet(ncolor)))
colorbar('Location','North','XTick',log10([2.5 5 10 20 50 100]),'XTickLabel',[2.5 5 10 20 50 100])
set(gca,'Visible','off')