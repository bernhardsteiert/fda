% Figure 1b (eyecatcher): Traces of BTC 100ng/ml in 184A1
addpath('./Functions/')

isite = 64;
% isite = 4;

load(['./Workspaces/site_' num2str(isite)])
c_signal = log10(intensity);

% % Interpolate (looks a bit weird --> switched off)
% c_signal(isinf(c_signal)) = nan;
% for i = 1:size(c_signal,2)
%     if sum(isnan(c_signal(:,i))) > length(c_signal(:,i))-2
%         c_signal(:,i) = 0;
%     end
%     c_signal(:,i) = interp1(timestamp(~isnan(c_signal(:,i))),c_signal(~isnan(c_signal(:,i)),i),timestamp);
%     vec = ~isnan(c_signal(:,i))';
%     rl = find(vec ~= [vec(2:end), vec(end)+1]);
%     data =  vec(rl);
%     rl(2:end) = rl(2:end) - rl(1:end-1);
%     if ~data(1)
%         c_signal(1:rl(1),i) = c_signal(rl(1)+1,i);
%     end
%     if ~data(end)
%         c_signal(end-rl(end)+1:end,i) = c_signal(end-rl(end),i);
%     end
% end

time_range = [50 1000];

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

figure;

s = siteprop(isite);

first_n = 7; % Plot first_n traces colored (2 times)

tmp_puls_strengths = edge_snr_score_pw_distdur(isite,[],0,1/120,'harm_basis_130722_corrected_retracked_all_cleaned_late',1);
[tmp ind_tmp_pul_str] = sort(tmp_puls_strengths);
ind_isite = [ind_tmp_pul_str(end:-1:end-first_n+1) ind_tmp_pul_str(round(linspace(1,(length(ind_tmp_pul_str)-first_n),first_n)))];

plot(repmat(timestamp(range_ind),1,size(c_signal,2)),c_signal(range_ind,:),'g','color',[0.7 0.7 0.7])
hold on
plot(repmat(timestamp(range_ind),1,2*first_n),c_signal(range_ind,ind_isite))
plot(timestamp(range_ind),c_signal(range_ind,ind_isite(1)),'color','k','LineWidth',2)
title([s.lig_name num2str(s.lig_dose) ' ng/ml'])

xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');

ylim = [-1 1]*.04;
xlim = get(gca,'XLim');
plot([200 200],ylim,'k--')
text(125,.033,{'non-stationary','(deterministic)'},'HorizontalAlignment','center')
text((xlim(2)-150)/2+150,.033,{'stationary','(stochastic)'},'HorizontalAlignment','center')

set(gca,'XLim',time_range,'YLim',ylim)
set(gca,'XTick',100:100:xlim(2))