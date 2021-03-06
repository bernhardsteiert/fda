% Figure 2b: Traces of all ligands at 100ng/ml in 184A1
addpath('./Functions/')

% sites_all = [17 37 44 4 57 64];
sites_all = [17 37 44 64 57 4];

times = cell(0);
signals = cell(0);
celltype = [];
legstr = {};

for isite = sites_all
    load(['./Workspaces/site_' num2str(isite)])
    times{end+1} = timestamp;  
    signals{end+1} = log10(intensity);
    celltype = [celltype ones(1,size(intensity,2))*isite];
    s = siteprop(isite);
    legstr{end+1} = s.lig_name;
end

timestamp = times{1}; % same time sampling for all data sets
c_signal = cell2mat(signals);

time_range = [50 605];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max+1;

figure

ncols = 3;
nrows = 2;

colmap = [linspace(0,1,length(sites_all)+1)' ones(length(sites_all)+1,1) ones(length(sites_all)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));

for iplot = 1:length(sites_all)-1
    subplot(nrows,ncols,iplot+1)
%     rectangle('Position',[65 0.02 50 0.015],'Curvature',[.2 .2],'FaceColor',colmap(iplot,:))
    hold on
    
    isite = sites_all(iplot);
    s = siteprop(isite);

    first_n = 7; % Plot first_n traces colored (2 times)

    tmp_puls_strengths = edge_snr_score_pw_distdur(isite,[],0,1/120,'harm_basis_130722_corrected_retracked_all_cleaned_late',1);
    [tmp ind_tmp_pul_str] = sort(tmp_puls_strengths);
    ind_isite = [ind_tmp_pul_str(end:-1:end-first_n+1) ind_tmp_pul_str(round(linspace(1,(length(ind_tmp_pul_str)-first_n),first_n)))];

    c_signal_single = c_signal(:,celltype == isite);
    plot(repmat(timestamp(range_ind),1,size(c_signal_single,2)),c_signal_single(range_ind,:),'g','color',[0.7 0.7 0.7])
    plot(repmat(timestamp(range_ind),1,2*first_n),c_signal_single(range_ind,ind_isite))
    title(s.lig_name)

    if iplot == 5
        xlabel('time [min]')
    end
    if iplot == 1
        ylabel('log_{10} FOXO3a [Cyt/Nuc]')
    end

    ylim = [-1 1]*.05;
%     plot([200 200],ylim,'k--')
%     text(125,.033,{'non-stationary','(deterministic)'},'HorizontalAlignment','center')
%     text(355,.033,{'stationary','(stochastic)'},'HorizontalAlignment','center')

    set(gca,'XLim',[50 600],'YLim',ylim);
%     if mod(iplot,2)==0
%         set(gca,'YAxisLocation','right');
%     end
    set(gca,'XTick',120:100:520,'XTickLabel',0:100:400,'YTick',[-0.04 -0.02 0 0.02 0.04])
end

subplot(nrows,ncols,1)
hold on
% resort = [1 2 3 4 6 5];
resort = [1 2 3 6 4 5];
legh = [];
for iplot = 1:length(sites_all)
    isite = sites_all(resort(iplot));
    legh = [legh plot(timestamp(range_ind),nanmean(c_signal(range_ind,celltype == isite),2),'Color',colmap(iplot,:))];
end
title('Averaged')
set(gca,'XLim',time_range,'YLim',ylim/2)
%set(gca,'XTick',0:200:time_range(2),'YTick',[-0.01 0 0.01])
set(gca,'XTick',120:100:520,'XTickLabel',0:100:400,'YTick',[-0.04 -0.02 0 0.02 0.04])
legend(legh,legstr{resort})