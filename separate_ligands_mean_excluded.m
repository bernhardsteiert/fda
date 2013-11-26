%% Look at qualitatively different time-courses and use fPCA to reseparate them
% Standard fPCA: mean + sum(coef*harm)
% Problem: Mean behavior is meaningless when mean is over different stimuli;
% Solution: Do fPCA without for mean behavior first ...

close all
clear all
clc

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
addpath(grabdataPath)

% Possible stages/sites:
% sites = [4 17 37 44 57 64];
% input_names = {'EGF','IGF','HRG','HGF','EPR','BTC'};
sites = [17 57 64];
input_names = {'IGF','EPR','BTC'};
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

%% Generate spline fits to individual data sets with nbasis basis functions
close all

nbasis = 50;
% time_range = [min(timestamp) max(timestamp)];
time_range = [50 650];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
smoothed_data = smooth_basis(timestamp(range_ind),c_signal(range_ind,:),basis);

f = figure;
set(f,'DefaultAxesColorOrder',jet(size(c_signal,2)))
hold on

plot(smoothed_data)
plot(timestamp(range_ind),c_signal(range_ind,:),'o')

%% Make FPCA with data generated in previous block - wip
close all

nharm = 8;
c_signal_pcastr = pca_fd(smoothed_data, nharm);

plot_pca_fd(c_signal_pcastr, 1, 0)

%% Plot: %variance explained vs. #basis functions
close all

thres_var = 0.9;

cumprobs = cumsum([0;c_signal_pcastr.varprop]);
[tmp thres_ind] = min(abs(cumprobs - thres_var));

plot(0:length(c_signal_pcastr.varprop),cumprobs)
hold on
plot(0:thres_ind-1,ones(1,thres_ind)*thres_var,'--')
plot([thres_ind-1 thres_ind-1],[0 thres_var],'--')

xlabel('fPCA basis functions')
ylabel('cumulative variance explained')

fprintf('To explain %s variance, use %i fPCA basis functions.\n\n',num2str(thres_var,3),thres_ind-1);

%% Plot: Eigenvalues linear trend
close all

linearrange = 8:15;

axb = polyfit(linearrange,log10(c_signal_pcastr.values(linearrange))',1);
plot(1:max(linearrange),log10(c_signal_pcastr.values(1:max(linearrange))),'ko-')
hold on
plot(1:max(linearrange),(1:max(linearrange))*axb(1) + axb(2),'k--')

xlabel('Eigenvalue Number')
ylabel('log10(Eigenvalue)')

%% Plot: Eigenfunctions
% close all
figure

rowstocols = 0.5;
nrows = ceil((nharm+1)^rowstocols);
ncols = ceil((nharm+1) / nrows);

time_range = [50 650];

flipharm = ones(1,nharm);
flipharm(1:8) = [1 -1 1 -1 -1 1 1 1];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;
times_fine = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),501);

harm_eval = 2 * repmat(sqrt(c_signal_pcastr.values(1:nharm))'.*flipharm(1:nharm),length(times_fine),1) .* eval_fd(c_signal_pcastr.harmfd,times_fine);

subplot(nrows,ncols,1)
plot(times_fine,eval_fd(c_signal_pcastr.meanfd,times_fine))
xlabel('Mean')
set(gca,'XLim',time_range)
set(gca,'YLim',[min(min(harm_eval)) max(max(harm_eval))])
hold on
plot(time_range,[0 0],'--')

for iplot = 1:nharm
    subplot(nrows,ncols,iplot+1)
    
    plot(times_fine,harm_eval(:,iplot))
    xlabel(['Harmonic ' num2str(iplot)])
    set(gca,'XLim',time_range)
    set(gca,'YLim',[min(min(harm_eval)) max(max(harm_eval))])
    
    hold on
    plot(time_range,[0 0],'--')
end

%% Plot: Harmonic scores PCi vs. PCj
close all

figure
color = lines(length(signals));
hold on

unitypes = unique(celltype);
linewidth = 1;
markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};

pcs = [1 2];
plotfat = [17 57];

for ilig = 1:length(signals)
    plot(c_signal_pcastr.harmscr(celltype == unitypes(ilig),pcs(1)),c_signal_pcastr.harmscr(celltype == unitypes(ilig),pcs(2)),markers{ilig},'Color',color(ilig,:),'LineWidth',linewidth+sum(unitypes(ilig)==plotfat))
end

xlabel(['PC ' num2str(pcs(1))])
ylabel(['PC ' num2str(pcs(2))])
legend(input_names)

%% Plot: Triagonal Matrix of PCs
close all

max_pc = 8;
% plotfat = [17 57];
plotfat = [];

flipharm = ones(1,nharm);
flipharm(1:8) = [1 -1 1 -1 -1 1 1 1];

harm_scores_flipped = c_signal_pcastr.harmscr(:,1:max_pc) .* repmat(flipharm(1:max_pc),size(c_signal_pcastr.harmscr(:,max_pc),1),1);

figure
color = lines(length(signals));
legendstyles = nan(1,length(signals));
hold on

unitypes = unique(celltype);
linewidth = 1;
markers = {'+','o','*','x','s','d','^','v','>','<','p','h','.'};

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
            legendstyles(ilig) = plot(harm_scores_flipped(celltype == unitypes(ilig),icol),harm_scores_flipped(celltype == unitypes(ilig),irow+1),markers{ilig},'Color',color(ilig,:),'LineWidth',linewidth+sum(unitypes(ilig)==plotfat));
        end
        
        xrange = max(harm_scores_flipped(:,icol)) - min(harm_scores_flipped(:,icol));
        yrange = max(harm_scores_flipped(:,irow+1)) - min(harm_scores_flipped(:,irow+1));
        scalefac = .05;
        
        set(gca,'XLim',[min(harm_scores_flipped(:,icol)) - xrange*scalefac max(harm_scores_flipped(:,icol)) + xrange*scalefac])
        set(gca,'YLim',[min(harm_scores_flipped(:,irow+1)) - yrange*scalefac max(harm_scores_flipped(:,irow+1)) + yrange*scalefac])
        
    end
end

g = subplot(max_pc-1,max_pc-1,max_pc-1);
set(gca,'Visible','off')
legend(g,legendstyles,input_names)
