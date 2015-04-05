function gen_singlecelltraces(isite)
% Figure 1b (eyecatcher): Traces of BTC 100ng/ml in 184A1
addpath('./Functions/')
myextension = '130722_corrected_retracked_all_paper_cleaned';
% isite = 64;

% Change this path to match the location where you save ND file and would
% like to store PNG files
ndfilename = '130722.nd'; 
ndpathname = 'Z:\computation\Bernhard_Steiert\FOXO3a dynamics\Images and Data';
savepath  = 'C:\Users\SOMPONNAT\Dropbox (Somponnat workspace)\Somponnat paper\website\page 1 - response\PNGs for grids and popups\PNG_singlecelltraces';
prefix = ndfilename(1:(end-3));
[notp stagePos stageName channelnames] = readndfile(ndpathname,ndfilename);
tokens   = regexp(stageName{isite}, 'r(?<row>\d+)c(?<col>\d+)|r(?<row>\d+)_c(?<col>\d+)|R(?<row>\d+)C(?<col>\d+)|R(?<row>\d+)_C(?<col>\d+)','tokens');
if ~isempty(tokens)
    row = str2num(tokens{1}{1});
    col = str2num(tokens{1}{2});
else
    row = site;
    col = 1;
end

close all

load(['./Workspaces/site_' num2str(isite) '_' myextension])
c_signal = log10(intensity);

time_range = [50 605];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

figure;

s = siteprop(isite);

first_n = 7; % Plot first_n traces colored (2 times)

tmp_puls_strengths = edge_snr_score_pw_distdur(isite,myextension,0,1/120,'harm_basis_130722_corrected_retracked_all_cleaned_late',1);
[tmp ind_tmp_pul_str] = sort(tmp_puls_strengths);
ind_isite = [ind_tmp_pul_str(end:-1:end-first_n+1) ind_tmp_pul_str(round(linspace(1,(length(ind_tmp_pul_str)-first_n),first_n)))];

plot(repmat(timestamp(range_ind),1,size(c_signal,2)),c_signal(range_ind,:),'g','color',[0.7 0.7 0.7])
hold on
plot(repmat(timestamp(range_ind),1,2*first_n),c_signal(range_ind,ind_isite))
% plot(timestamp(range_ind),nanmean(c_signal(range_ind,:),2),'color','k','LineWidth',2)

xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');

ylim = [-1 1]*.04;
plot([200 200],ylim,'k--')
text(125,.033,{'non-stationary','(deterministic)'},'HorizontalAlignment','center')
text(400,.033,{'stationary','(stochastic)'},'HorizontalAlignment','center')

set(gca,'XLim',[50 600],'YLim',ylim)
set(gca,'XTick',120:100:520,'XTickLabel',0:100:400)

set(gcf,'Position',[100 100 700 400]);
%set(gca,'Position',[0.1 0.1 0.95 0.95],'Color','w','Visible','off');
export_fig(fullfile(savepath,['s_r' num2str(row) 'c' num2str(col) '.png']),'-nocrop','-png','-transparent');
