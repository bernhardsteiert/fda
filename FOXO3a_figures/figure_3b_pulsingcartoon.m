% Figure 3B: Explanation of pulsing score calculation
addpath('./Functions/')

close all

load('./Workspaces/harm_basis_fPCA_5basis_noFGF_newBTC_rot')
harm_basis_fPCA = harm_basis;
load('./Workspaces/harm_basis_130722_corrected_retracked_all_cleaned_late_newBTC')
load('./Workspaces/scores_early_5basis_noFGF_newBTC')

myextension = '130722_corrected_retracked_all_cleaned';

sites_all = [17 57 64];
% sigs = [26 6 17]; % IGF 26 looks OK; EPR 6 looks OK; BTC 32 has transient shape + pulsing (75 is maybe better; 17 has some of everything)
sigs = [26 17 99];
colind = [1 6 5];

colmap = [linspace(0,1,7)' ones(7,1) ones(7,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
markers = {'o','s','v','d','^','>'};
linewidth = 2;

xfac = 1;
yfac = 1;
fontsize = 24;

legstr = cell(1,length(sites_all));
c_signal_single = [];
scores_single = nan(5,length(sigs));

for icount = 1:length(sites_all)
    isite = sites_all(icount);
    load(['./Workspaces/site_' num2str(isite) '_' myextension])
    c_signal_single(:,icount) = log10(intensity(:,sigs(icount)));
    scores_tmp = scores_early(:,celltypes==isite);
    scores_single(:,icount) = scores_tmp(:,sigs(icount));
    s = siteprop(isite);
    legstr{icount} = s.lig_name;
end

c_signal_single(isinf(c_signal_single)) = nan;
for i = 1:size(c_signal_single,2)
    if sum(isnan(c_signal_single(:,i))) > length(c_signal_single(:,i))-2
        c_signal_single(:,i) = 0;
    end
    c_signal_single(:,i) = interp1(timestamp(~isnan(c_signal_single(:,i))),c_signal_single(~isnan(c_signal_single(:,i)),i),timestamp);
    vec = ~isnan(c_signal_single(:,i))';
    rl = find(vec ~= [vec(2:end), vec(end)+1]);
    data =  vec(rl);
    rl(2:end) = rl(2:end) - rl(1:end-1);
    if ~data(1)
        c_signal_single(1:rl(1),i) = c_signal_single(rl(1)+1,i);
    end
    if ~data(end)
        c_signal_single(end-rl(end)+1:end,i) = c_signal_single(end-rl(end),i);
    end
end

nbasis = 20;
time_range = [50.7 197.8];
[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;
times_fine = linspace(time_range(1),time_range(2),501);
basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
smoothed_data = smooth_basis(timestamp(range_ind),c_signal_single(range_ind,:),basis);
smoothed_eval = eval_fd(smoothed_data,times_fine);
times_fine(end+1) = time_range(2);
smoothed_eval(end+1,:) = smoothed_eval(end,:);

figure
setFigure(gcf,xfac,yfac,fontsize)

hold on
legh = [];
for iplot = 1:size(c_signal_single,2)
    legh = [legh plot(timestamp,c_signal_single(:,iplot),markers{colind(iplot)},'Color',colmap(colind(iplot),:),'LineWidth',linewidth)];
    plot(timestamp,timestamp*0,'k:','LineWidth',linewidth)
    plot(times_fine,smoothed_eval(:,iplot),'Color',colmap(colind(iplot),:),'LineWidth',linewidth)
end
set(gca,'XLim',[50 200],'XTick',60:30:180,'XTickLabel',-60:30:60)

h = legend(legh,legstr,'Location','NorthWest');
ch = get(h,'child');
for ileg = 1:length(ch)/3
    ilegch = (ileg-1)*3+2;
    set(ch(ilegch),'LineStyle','-','Color',colmap(colind(end-ileg+1),:)); 
end

xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');

figure
setFigure(gcf,xfac,yfac,fontsize)

hold on

legh = [];
for isite = 1:size(c_signal_single,2)
    legh = [legh plot(timestamp,c_signal_single(:,isite),[markers{colind(isite)} '-'],'Color',colmap(colind(isite),:),'LineWidth',linewidth)];
    plot(timestamp,timestamp*0,'k:','LineWidth',linewidth)
end
set(gca,'XLim',[50 200],'XTick',60:30:180,'XTickLabel',-60:30:60)
legend(legh,legstr,'Location','NorthWest')

xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');

figure

nrows = 3;
ncols = 2;

hold on

legh = [];
for isite = 1:size(c_signal_single,2)
    legh = [legh plot(timestamp,c_signal_single(:,isite),[markers{colind(isite)} '-'],'Color',colmap(colind(isite),:))];
    plot(timestamp,timestamp*0,'k:')
end

legend(legh,legstr)

set(gca,'XLim',[50 1000],'XTick',120:100:920,'XTickLabel',0:100:800)
% plot([200 200],get(gca,'YLim'),'k--')

xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');

time_range = [50.7 197.8];
times_fine = linspace(time_range(1),time_range(2),501);

basis_eval = eval_basis(harm_basis_fPCA,times_fine);
times_fine(end+1) = 200;
basis_eval(end+1,:) = basis_eval(end,:);

axpos = get(gca,'Position');
xlim = get(gca,'XLim');
annotation('rectangle',[axpos(1:2) (times_fine(end)-xlim(1))./range(xlim) * axpos(3) axpos(4)])


figure
setFigure(gcf,xfac,yfac,fontsize)

nrows = 3;
ncols = 2;

hold on

legh = [];
for isite = 1:size(c_signal_single,2)
    legh = [legh plot(timestamp,c_signal_single(:,isite),[markers{colind(isite)} '-'],'Color',colmap(colind(isite),:),'LineWidth',linewidth)];
    plot(timestamp,timestamp*0,'k:','LineWidth',linewidth)
end

legend(legh,legstr,'Location','SouthEast')

set(gca,'XLim',[200 1000],'XTick',220:100:920,'XTickLabel',100:100:800)
% plot([200 200],get(gca,'YLim'),'k--')

xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');

time_range = getbasisrange(harm_basis);
[tmp range_ind_min] = min(abs(timestamp - time_range(1) - 5));
[tmp range_ind_max] = min(abs(timestamp - time_range(2) + 5));
range_ind = range_ind_min:range_ind_max;
times_fine_late = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),201);
smoothed_additional = smooth_basis(timestamp(range_ind),c_signal_single(range_ind,:),harm_basis);
harm_eval = eval_basis(harm_basis,timestamp(range_ind));
harm_eval_fine = eval_basis(harm_basis,times_fine_late);
fitcoef = getcoef(smoothed_additional);
data_fpca_repr = fitcoef'*harm_eval';
data_fpca_repr_fine = fitcoef'*harm_eval_fine';
c_signal_woNharm = c_signal_single(range_ind,:)-data_fpca_repr';

figure
setFigure(gcf,xfac,yfac,fontsize)

hold on
legh = [];
for iplot = 1:size(c_signal_single,2)
    legh = [legh plot(timestamp,c_signal_single(:,iplot),markers{colind(iplot)},'Color',colmap(colind(iplot),:),'LineWidth',linewidth)];
    plot(timestamp,timestamp*0,'k:','LineWidth',linewidth)
    plot(times_fine_late,data_fpca_repr_fine(iplot,:),'Color',colmap(colind(iplot),:),'LineWidth',linewidth)
end
set(gca,'XLim',[200 1000],'XTick',220:100:920,'XTickLabel',100:100:800)

h = legend(legh,legstr,'Location','SouthEast');
ch = get(h,'child');
for ileg = 1:length(ch)/3
    ilegch = (ileg-1)*3+2;
    set(ch(ilegch),'LineStyle','-','Color',colmap(colind(end-ileg+1),:)); 
end

xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');

figure
setFigure(gcf,xfac,yfac,fontsize)

nrows = 3;
ncols = 2;

hold on

legh = [];
for isite = 1:size(c_signal_single,2)
    legh = [legh plot(timestamp(range_ind),c_signal_woNharm(:,isite),[markers{colind(isite)} '-'],'Color',colmap(colind(isite),:),'LineWidth',linewidth)];
    plot(timestamp,timestamp*0,'k:','LineWidth',linewidth)
end

legend(legh,legstr,'Location','SouthEast')

set(gca,'XLim',[200 1000],'XTick',220:100:920,'XTickLabel',100:100:800)
set(gca,'YLim',[-.015 .015])
% plot([200 200],get(gca,'YLim'),'k--')

xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');

figure
setFigure(gcf,xfac,yfac,fontsize)

hold on
legh = [];
for iplot = 1:size(c_signal_single,2)
    legh = [legh plot(timestamp,c_signal_single(:,iplot),markers{colind(iplot)},'Color',colmap(colind(iplot),:),'LineWidth',linewidth)];
    plot(timestamp,timestamp*0,'k:','LineWidth',linewidth)
    plot(times_fine,basis_eval*scores_single(:,iplot),'Color',colmap(colind(iplot),:),'LineWidth',linewidth)
end
set(gca,'XLim',[50 200],'XTick',60:30:180,'XTickLabel',-60:30:60)

h = legend(legh,legstr,'Location','NorthWest');
ch = get(h,'child');
for ileg = 1:length(ch)/3
    ilegch = (ileg-1)*3+2;
    set(ch(ilegch),'LineStyle','-','Color',colmap(colind(end-ileg+1),:)); 
end

xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');

figure
setFigure(gcf,xfac*1.5,yfac,fontsize)

nbasis = round(length(range_ind)/3.15);
basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
smoothed_data_woNharm = smooth_basis(timestamp(range_ind),c_signal_woNharm,basis);
c_smoothed_eval = eval_fd(smoothed_data_woNharm,linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),201));
xs = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),201);

noise_thres = .6;
range_smoothed_in = 1/120;

c_smoothed_eval_data = eval_fd(smoothed_data_woNharm,timestamp(range_ind));
rss_spline = nansum((c_signal_woNharm - c_smoothed_eval_data).^2,1);

nEdges = [];
SNR = [];
amp = [];

peakdur_mean = [];
peakdur_std = [];
peakdis_mean = [];
peakdis_std = [];

peakdur = [];
peakdis = [];

for isig = 1:size(c_signal_woNharm,2)
    [pks,locs] = findpeaks(c_smoothed_eval(:,isig));
    [pks2,locs2] = findpeaks(-c_smoothed_eval(:,isig));

    all_locs = [locs; locs2];
    all_pks = [pks; -pks2];
    all_type = [ones(1,length(pks)) -ones(1,length(pks2))];

    [locs_sorted ind_locs_sorted] = sort(all_locs);
    % --> all_type(ind_locs_sorted) is alternating by construction
    locs_sorted = [1; locs_sorted; length(xs)];
    types_sorted = all_type(ind_locs_sorted);
    if ~isempty(types_sorted)
        types_sorted = [-types_sorted(1) types_sorted -types_sorted(end)];
    else
        types_sorted = nan;
    end

%         range_smoothed = max(c_smoothed_eval(:,isig))-min(c_smoothed_eval(:,isig)); % relative
%         range_smoothed = 1/60; % absolute --> edge only counted when > 0.01
%         range_smoothed = 1/120; % absolute --> edge only counted when > 0.005

    % This is new compared to before:
%         range_smoothed_in = max([(max(c_smoothed_eval(:,isig))-min(c_smoothed_eval(:,isig)))/2 range_smoothed_in]);

    candidate_left = [c_smoothed_eval(1,isig); all_pks(ind_locs_sorted(1:end))];
    candidate_right = [all_pks(ind_locs_sorted(1:end)); c_smoothed_eval(end,isig); nan];

    icl = 1;
    ind_edge_start = [];
    ind_edge_end = [];
    while icl < length(candidate_right)
        if icl < length(candidate_right)-1
        testcan = types_sorted(icl)*(candidate_left(icl)-candidate_right([icl; icl+2]));

            if ~(icl == 1 && testcan(1) < noise_thres*range_smoothed_in && testcan(1) < testcan(2))

                if testcan(1) < testcan(2) && types_sorted(icl+1)*(candidate_left(icl+1)-candidate_right(icl+1)) < noise_thres*range_smoothed_in/3
                    % Three edges connected
                    ind_edge_start = [ind_edge_start icl];
                    ind_edge_end = [ind_edge_end icl+2];
                    icl = icl + 3;
                else
                    % Usual case; Neighboring peaks are connected by one edge
                    ind_edge_start = [ind_edge_start icl];
                    ind_edge_end = [ind_edge_end icl];
                    icl = icl+1;
                end
            else
                icl = icl + 1;
            end

        elseif types_sorted(icl)*(candidate_left(icl)-candidate_right(icl)) > noise_thres*range_smoothed_in
            % Usual case; Neighboring peaks are connected by one edge
            ind_edge_start = [ind_edge_start icl];
            ind_edge_end = [ind_edge_end icl];
            icl = icl+1;
        else
            icl = icl + 1;
        end
    end

    edge_heights = candidate_left(ind_edge_start)-candidate_right(ind_edge_end);
    ispeak = abs(edge_heights) > noise_thres*range_smoothed_in;

    peak_duration = [];
    peak_distance = [];
    ind_peaks = find(ispeak)';
    
    mycount1 = 1;
    mycount2 = 1;
    for i = ind_peaks
        if i < length(ispeak)
            if ispeak(i+1)
                peak_duration = [peak_duration xs(locs_sorted(ind_edge_end(i+1)+1))-xs(locs_sorted(ind_edge_start(i)))];
                if isig == 3 && xs(locs_sorted(ind_edge_end(i+1)+1)) < 1000
                    subplot(1,3,[1 2])
                    mypos = get(gca,'Position');
                    annotation('doublearrow',(([xs(locs_sorted(ind_edge_start(i))) xs(locs_sorted(ind_edge_end(i+1)+1))]-200)./800)*mypos(3)+mypos(1),(mypos(2)*1.1+.05*mod(mycount1,2)*mypos(4))*[1 1])
                    mycount1 = mycount1 + 1;
                end
            end

            tmpind = find(ind_peaks > i);
            tmpind2 = find(types_sorted(ind_edge_start(ind_peaks(tmpind))) == types_sorted(ind_edge_start(i)));
            if ~isempty(tmpind2)
                peak_distance = [peak_distance xs(locs_sorted(ind_edge_start(ind_peaks(tmpind(tmpind2(1))))))-xs(locs_sorted(ind_edge_start(i)))];
                if isig == 3 && xs(locs_sorted(ind_edge_start(ind_peaks(tmpind(tmpind2(1)))))) < 1000
                    subplot(1,3,[1 2])
                    mypos = get(gca,'Position');
                    annotation('doublearrow',(([xs(locs_sorted(ind_edge_start(i))) xs(locs_sorted(ind_edge_start(ind_peaks(tmpind(tmpind2(1))))))]-200)./800)*mypos(3)+mypos(1),(mypos(4)+mypos(2)*.9-.05*mod(mycount2,2)*mypos(4))*[1 1])
                    mycount2 = mycount2 + 1;
                end
            end

        end
    end

    peakdur = [peakdur peak_duration];
    peakdis = [peakdis peak_distance];

    peakdur_std = [peakdur_std std(peak_duration)];
    peakdis_std = [peakdis_std std(peak_distance)];

    tmp = mean(peak_duration); peakdur_mean = [peakdur_mean min([tmp 300])];
    tmp = mean(peak_distance); peakdis_mean = [peakdis_mean min([tmp 300])];

%         tmp = min(peak_duration); peakdur_mean = [peakdur_mean min([tmp 300])]; %% !!!
%         tmp = min(peak_distance); peakdis_mean = [peakdis_mean min([tmp 300])]; %% !!!

    nEdges = [nEdges sum(ispeak)];
    range_smoothed = max(c_smoothed_eval(:,isig))-min(c_smoothed_eval(:,isig));
    SNR = [SNR range_smoothed./sqrt(rss_spline(isig)./(size(c_signal_woNharm,1)-1))];
%         if ~isnan(median(edge_heights(ispeak)))
%             amp = [amp median(abs(edge_heights(ispeak)))];
%         else
        amp = [amp range_smoothed]; % Only amplitude
%         end

    if isig == 3
%                 close all
%                 figure
        subplot(1,3,[1 2])
        plot(xs,c_smoothed_eval(:,isig),'k--')
        hold on
        plot(timestamp(range_ind),c_signal_woNharm(:,isig),'^')
        hold on
        mycolor = lines(length(ispeak));
        for i = 1:length(ind_edge_start)
            myrange = find(candidate_left(ind_edge_start(i))==c_smoothed_eval(:,isig)):find(candidate_right(ind_edge_end(i))==c_smoothed_eval(:,isig));
            if ispeak(i)
                plot(xs(myrange),c_smoothed_eval(myrange,isig),'Color',mycolor(i,:),'LineWidth',2)
            end
        end
        plot(get(gca,'XLim'),[max(c_smoothed_eval(:,isig)) max(c_smoothed_eval(:,isig))],'k:')
        plot(get(gca,'XLim'),[min(c_smoothed_eval(:,isig)) min(c_smoothed_eval(:,isig))],'k:')
        set(gca,'XLim',[200 1000],'XTick',220:100:920,'XTickLabel',100:100:800)
        subplot(1,3,3)
        set(gca,'Visible','off')
        text(.1,.9,sprintf('nEdges = %g',nEdges(end)))
        text(.1,.8,sprintf('SNR = %g',SNR(end)))
        text(.1,.7,sprintf('Amplitude = %g',amp(end)))
        text(.1,.6,sprintf('Peak dur mean = %g',peakdur_mean(end)))
        text(.1,.5,sprintf('Peak dur std = %g',peakdur_std(end)))
        text(.1,.4,sprintf('Peak dis mean = %g',peakdis_mean(end)))
        text(.1,.3,sprintf('Peak dis std = %g',peakdis_std(end)))
        1;

%             tmp = [];
%             for i = 1:length(ind_edge_start)
%                 if ispeak(i)
%                     tmp = padconcatenation(tmp,find(candidate_left(ind_edge_start(i))==c_smoothed_eval(:,isig)):find(candidate_right(ind_edge_end(i))==c_smoothed_eval(:,isig)),1);
%                 end
%             end
    end

end
