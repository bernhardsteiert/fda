%% Data-sets with different drug and ligand concentrations

close all
clear all
clc

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)

grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
addpath(grabdataPath)

dataPath = '2D_dose_response_drugsVSEGF_130903';

% site 20 (row 2 - col 1) is missing ...
sites = [1:19 21:70];
sites_for_harmonics = [1:19 21:70];
lig_name = 'EGF';


times = cell(0);
signals = cell(0);
celltype = [];

for isite = sites
    if exist(remotepath,'dir')
        [times{end+1},intensity] = grabdata(isite,dataPath);
    else
        load(['./Workspaces/site_' num2str(isite) '_only_' lig_name])
        times{end+1} = timestamp;
    end

    log_trafo = 0; % log-transform signal

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

%% Reproduce plot from Bernhard Kraemer
close all

site_lig_dose = [];
site_inh_name = {};
site_inh_dose = [];

for isite = sites_for_harmonics   % only plot ligands used for harmonics
    s = siteprop_drug(isite);
    site_lig_dose = [site_lig_dose s.lig_dose];
    site_inh_name{end+1} = s.inh_name;
    site_inh_dose = [site_inh_dose s.inh_dose];
end

figure

nrows = 7;
ncols = 10;

for iplot = sites
    isite = iplot;
    % Switch every second row
    inh_index = ceil(iplot/10);
    dose_index = mod(iplot-1,10)+1;
    if ~mod(inh_index,2)
        dose_index = 11-dose_index;
    end
    isite = 10*(inh_index-1)+dose_index;
    
    subplot(nrows,ncols,isite)
    hold on
    
    plot(repmat(timestamp,1,sum(celltype == iplot)),c_signal(:,celltype == iplot),'g','color',[0.7 0.7 0.7])
    plot(timestamp,nanmean(c_signal(:,celltype == iplot),2),'color','k','LineWidth',2)
    prop_ind = (iplot == sites);
    title([lig_name ' ' num2str(site_lig_dose(prop_ind)) site_inh_name{prop_ind} num2str(site_inh_dose(prop_ind))])
    
    ylim = [-1 1]*.04;
    if ~log_trafo
        ylim = 10.^ylim;
    end
    set(gca,'XLim',[0 500],'YLim',ylim)

end
