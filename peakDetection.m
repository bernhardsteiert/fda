%% Examples: isig = 5 --> should have low score; isig = [1 2 3] should have high scores
isig = 5;

load('harm_basis.mat') % Contains only harm_basis from all data-sets

time_range = getbasisrange(harm_basis);
[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;
if timestamp(range_ind(1)) < time_range(1)
    range_ind = range_ind(2:end);
end
if timestamp(range_ind(end)) > time_range(2)
    range_ind = range_ind(1:end-1);
end

smoothed_additional = smooth_basis(timestamp(range_ind),c_signal_single(range_ind,:),harm_basis);

harm_eval = eval_basis(harm_basis,timestamp(range_ind));

fitcoef = getcoef(smoothed_additional);
data_fpca_repr = fitcoef'*harm_eval';

c_signal_woNharm = c_signal_single(range_ind,:)-data_fpca_repr';


MiddleToTop = 1; % 1 = assugb the middle cluster to top group, 0 = assign the middle cluster to bottom group
showPlots = 1; % change to 1 if needing to visualize the peak detection
midHgating = -0.01; % x the median of lowest peak cluster
delayGate = 0.5; % fraction of height that must decay to consider as peak tail

divisionTime = 10;

bsfindTruePeaks(timestamp(range_ind),c_signal_woNharm(:,end-isig+1),showPlots,MiddleToTop,midHgating,delayGate,[timestamp(range_ind(1)) timestamp(range_ind(end))],divisionTime);


%%
isig =  2;
close all

noise_thres = .4;

nbasis = 20;
basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
smoothed_data_woNharm = smooth_basis(timestamp(range_ind),c_signal_woNharm,basis);
c_smoothed_eval = eval_fd(smoothed_data_woNharm,linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),201));
c_smoothed_eval_data = eval_fd(smoothed_data_woNharm,timestamp(range_ind));

rss_spline = nansum((c_signal_woNharm - c_smoothed_eval_data).^2,1);

xs = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),201);
plot(xs,c_smoothed_eval(:,end-isig+1))
hold on
plot(timestamp(range_ind),c_signal_woNharm(:,end-isig+1),'-x')

[pks,locs] = findpeaks(c_smoothed_eval(:,end-isig+1));
% plot(xs(locs),pks,'gv')
[pks2,locs2] = findpeaks(-c_smoothed_eval(:,end-isig+1));
% plot(xs(locs2),-pks2,'r^')

all_locs = [locs; locs2];
all_pks = [pks; -pks2];
all_type = [ones(1,length(pks)) -ones(1,length(pks2))];

[locs_sorted ind_locs_sorted] = sort(all_locs);
% --> all_type(ind_locs_sorted) is alternating by construction
locs_sorted = [1; locs_sorted; size(c_signal,1)];
types_sorted = [-types_sorted(1) types_sorted -types_sorted(end)];

range_smoothed = max(c_smoothed_eval(:,end-isig+1))-min(c_smoothed_eval(:,end-isig+1));
candidate_left = [c_smoothed_eval(1,end-isig+1); all_pks(ind_locs_sorted(1:end))];
candidate_right = [all_pks(ind_locs_sorted(1:end)); c_smoothed_eval(end,end-isig+1)];
ispeak = abs(candidate_left-candidate_right) > noise_thres*range_smoothed;

blockstart = strfind([0 ~ispeak' 0], [0 1]);
blockend = strfind([0 ~ispeak' 0], [1 0]);
blocklength = blockend-blockstart;

% [maxlength blockpos] = max(blocklength);
bs = find(mod(blocklength,2));
while ~isempty(bs)
    b = bs(1);
    if blocklength(b) > 1
        candidate_left = candidate_left(setdiff(1:length(candidate_left),blockstart(b)+1:blockend(b)-1));
        candidate_right = candidate_right(setdiff(1:length(candidate_right),blockstart(b):blockend(b)-2));
    else
        % So not nicely done :(
        if blockstart(b) > 1
            if blockend(b) < length(candidate_left)+1
                candidate_left = candidate_left(setdiff(1:length(candidate_left),blockstart(b):blockend(b)));
                candidate_right = candidate_right(setdiff(1:length(candidate_right),blockstart(b)-1:blockend(b)-1));
            else
                candidate_left = candidate_left(setdiff(1:length(candidate_left),blockstart(b):blockend(b)));
                candidate_right = candidate_right(setdiff(1:length(candidate_right),blockstart(b):blockend(b)-1));
            end
        else
            candidate_left = candidate_left(setdiff(1:length(candidate_left),blockstart(b):blockend(b)-1));
            candidate_right = candidate_right(setdiff(1:length(candidate_right),blockstart(b)-1:blockend(b)-1));
        end
    end
    ispeak = abs(candidate_left-candidate_right) > noise_thres*range_smoothed;

    blockstart = strfind([0 ~ispeak' 0], [0 1]);
    blockend = strfind([0 ~ispeak' 0], [1 0]);
    blocklength = blockend-blockstart;
    bs = find(mod(blocklength,2));
end

nEdges = sum(ispeak)
SNR = range_smoothed./sqrt(rss_spline(end-isig+1)./size(c_signal_woNharm,1))

score = nEdges + SNR

options = statset('Display','off');
idx = kmeans(all_pks,3,'distance','sqEuclidean','emptyaction','drop','options',options,'Replicates',5);

matOrder(1,:) = [mean(all_pks(idx == 1)) mean(all_pks(idx == 2)) mean(all_pks(idx == 3))];
matOrder(2,:) = [std(all_pks(idx == 1)) std(all_pks(idx == 2)) std(all_pks(idx == 3))];
matOrder(3,:) = [1 2 3];
matOrder = sortrows(matOrder',1)';

TopP = all_locs(idx == matOrder(3,3));
MidP = all_locs(idx == matOrder(3,2));
BotP = all_locs(idx == matOrder(3,1));

plot(xs(TopP),all_pks(idx == matOrder(3,3)),'gv')
plot(xs(MidP),all_pks(idx == matOrder(3,2)),'ks')
plot(xs(BotP),all_pks(idx == matOrder(3,1)),'r^')

mycolor = lines(length(ispeak));
for i = find(ispeak)'
    myrange = find(candidate_left(i)==c_smoothed_eval(:,end-isig+1)):find(candidate_right(i)==c_smoothed_eval(:,end-isig+1));
    plot(xs(myrange),c_smoothed_eval(myrange,end-isig+1),'Color',mycolor(i,:),'LineWidth',2)
end