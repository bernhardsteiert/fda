%% Fit additional data-sets on harmonics generated from initial data-sets

close all
clear all
clc

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
addpath(grabdataPath)


sites = [17 41 42 44 57];
input_names = {'IGF','HGF-MEKi','HGF-AKTi','HGF','EPR'};
% sites_for_harmonics = [17 44 57];
sites_for_harmonics = [17 41 42 44 57];
% all ligands highest dose


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

nbasis = 40;
% time_range = [min(timestamp) max(timestamp)];
time_range = [50 650];

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

nharm = 12;
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

time_range = [50 650];

flipharm = ones(1,nharm);
% flipharm(1:8) = [1 -1 -1 -1 -1 1 1 1];

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

max_pc = 5;

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
