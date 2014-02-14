function res_out = res(p,yin)
% yin(:,1) = condition; yin(:,2) = time; yin(:,3:end) = data
% p(1:length(conds),:) = n-1 scaling parameters; p(length(conds):end) = ymean

conds = unique(yin(:,1));
times = unique(yin(:,2));
nobs = size(yin,2)-2;

res_out = [];
for ic = 1:length(conds)
    ycond = yin(yin(:,1) == conds(ic),:);
    for iobs = 1:nobs
        if ic > 1
            res_out = [res_out ycond(:,2+iobs) - p(length(conds):end,iobs) + p(ic-1,iobs)];
        else
            res_out = [res_out ycond(:,2+iobs) - p(length(conds):end,iobs)];
        end
    end
end

res_out = res_out(~isnan(res_out));