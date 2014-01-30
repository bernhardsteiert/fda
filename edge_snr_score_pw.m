function [radial_dist c_signal_woNharm range_ind nEdges SNR] = edge_snr_score_pw(isite,myextension,timeshift)
    if(~exist('myextension','var'))
        myextension = '';
    elseif(~isempty(myextension))
        myextension = ['_' myextension];
    end
    if(~exist('timeshift','var'))
        timeshift = 0;
    end
    
    load('harm_basis.mat') % Contains only harm_basis from all data-sets
    
    remotepath = mypath();
    
    warning('off','MATLAB:dispatcher:pathWarning')
    
    fdaMPath = [remotepath 'fda'];
    addpath(fdaMPath)

    grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
    addpath(grabdataPath)

    log_trafo = 1; % log-transform signal
    register = 1; % register IC50
    time_range = getbasisrange(harm_basis);
    
    if exist(remotepath,'dir')
        [timestamp,intensity] = grabdata_new(isite,myextension(2:end));
    else
        load(['./Workspaces/site_' num2str(isite) myextension])
    end

    if log_trafo
        c_signal = log10(intensity);
    else
        c_signal = intensity;
    end
    timestamp = timestamp - timeshift;
    
    if register
        c_signal = register_signal(c_signal,myextension(2:end));
    end
    c_signal = c_signal - repmat(nanmean(c_signal,2),1,size(c_signal,2));
    
    [tmp range_ind_min] = min(abs(timestamp - time_range(1)));
    [tmp range_ind_max] = min(abs(timestamp - time_range(2)));
    range_ind = range_ind_min:range_ind_max;
    
    if timestamp(range_ind(1)) < time_range(1)
        range_ind = range_ind(2:end);
    end
    if timestamp(range_ind(end)) > time_range(2)
        range_ind = range_ind(1:end-1);
    end
    
    smoothed_additional = smooth_basis(timestamp(range_ind),c_signal(range_ind,:),harm_basis);
    
    harm_eval = eval_basis(harm_basis,timestamp(range_ind));
    
    fitcoef = getcoef(smoothed_additional);
    data_fpca_repr = fitcoef'*harm_eval';
    
    c_signal_woNharm = c_signal(range_ind,:)-data_fpca_repr';
    
    % Generate spline fit to data-set given in isite (for remaining variation)
    nbasis = 20;
    basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
    smoothed_data_woNharm = smooth_basis(timestamp(range_ind),c_signal_woNharm,basis);
    c_smoothed_eval = eval_fd(smoothed_data_woNharm,linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),201));
    
    noise_thres = .6;
    
    c_smoothed_eval_data = eval_fd(smoothed_data_woNharm,timestamp(range_ind));
    rss_spline = nansum((c_signal_woNharm - c_smoothed_eval_data).^2,1);
    
    nEdges = [];
    SNR = [];
    amp = [];
    
    for isig = 1:size(c_signal_woNharm,2)
        [pks,locs] = findpeaks(c_smoothed_eval(:,isig));
        [pks2,locs2] = findpeaks(-c_smoothed_eval(:,isig));

        all_locs = [locs; locs2];
        all_pks = [pks; -pks2];
        all_type = [ones(1,length(pks)) -ones(1,length(pks2))];
        
        [locs_sorted ind_locs_sorted] = sort(all_locs);
        % --> all_type(ind_locs_sorted) is alternating by construction
        range_smoothed = max(c_smoothed_eval(:,isig))-min(c_smoothed_eval(:,isig));
        candidate_left = [c_smoothed_eval(1,isig); all_pks(ind_locs_sorted(1:end))];
        candidate_right = [all_pks(ind_locs_sorted(1:end)); c_smoothed_eval(end,isig)];
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

        nEdges = [nEdges sum(ispeak)];
        SNR = [SNR range_smoothed./sqrt(rss_spline(isig)./(size(c_signal_woNharm,1)-1))];
        amp = [amp range_smoothed]; % Only amplitude
        
    end
    
    max_nEdges = nbasis - 5;
    max_SNR = 40;
    max_amp = 0.05;
    max_pw = 0.1;
    if ~log_trafo
        max_amp = 10.^max_amp;
    end
    
    weight_edg = 2; % Weight between nEdges and SNR (Determine by PCA later?)
    weight_snr = 1;
    weight_amp = 1.5;
    weight_pw = 1.5;
    
    a = getcoef(smoothed_data_woNharm);
    b = a(2:end,:);
    
    pw = sqrt(sum((a(1:end-1,:)-b).^2,1));
    
%     princomp([log10(nEdges)' log10(SNR)' log10(amp)'])
    
    radial_dist = (nEdges/(max_nEdges*1.05)).^weight_edg .* (SNR/(max_SNR*1.05)).^weight_snr .* (pw/(max_pw*1.05)).^weight_pw;
    radial_dist = radial_dist.^.5; % To make distances smaller
end