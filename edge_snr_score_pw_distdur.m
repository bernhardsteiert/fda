% function [radial_dist c_signal_woNharm range_ind nEdges SNR amp pw peakdur peakdur_std peakdis peakdis_std] = edge_snr_score_pw_distdur(isite,myextension,timeshift)
function [radial_dist c_signal_woNharm range_ind nEdges SNR amp pw peakdur_mean peakdur_std peakdis_mean peakdis_std] = edge_snr_score_pw_distdur(isite,myextension,timeshift)
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
    xs = linspace(timestamp(range_ind(1)),timestamp(range_ind(end)),201);
    
    noise_thres = .6;
    
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
        types_sorted = [-types_sorted(1) types_sorted -types_sorted(end)];

%         range_smoothed = max(c_smoothed_eval(:,isig))-min(c_smoothed_eval(:,isig)); % relative
%         range_smoothed = 1/60; % absolute --> edge only counted when > 0.1
        range_smoothed = 1/120; % absolute --> edge only counted when > 0.05
        candidate_left = [c_smoothed_eval(1,isig); all_pks(ind_locs_sorted(1:end))];
        candidate_right = [all_pks(ind_locs_sorted(1:end)); c_smoothed_eval(end,isig); nan];

        icl = 1;
        ind_edge_start = [];
        ind_edge_end = [];
        while icl < length(candidate_right)
            if icl < length(candidate_right)-1
            testcan = types_sorted(icl)*(candidate_left(icl)-candidate_right([icl; icl+2]));

                if ~(icl == 1 && testcan(1) < noise_thres*range_smoothed && testcan(1) < testcan(2))

                    if testcan(1) < testcan(2) && types_sorted(icl+1)*(candidate_left(icl+1)-candidate_right(icl+1)) < noise_thres*range_smoothed/3
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

            elseif types_sorted(icl)*(candidate_left(icl)-candidate_right(icl)) > noise_thres*range_smoothed
                % Usual case; Neighboring peaks are connected by one edge
                ind_edge_start = [ind_edge_start icl];
                ind_edge_end = [ind_edge_end icl];
                icl = icl+1;
            else
                icl = icl + 1;
            end
        end

        ispeak = abs(candidate_left(ind_edge_start)-candidate_right(ind_edge_end)) > noise_thres*range_smoothed;

        peak_duration = [];
        peak_distance = [];
        ind_peaks = find(ispeak)';

        for i = ind_peaks
            if i < length(ispeak)
                if ispeak(i+1)
                    peak_duration = [peak_duration xs(locs_sorted(ind_edge_end(i+1)+1))-xs(locs_sorted(ind_edge_start(i)))];
                end

                tmpind = find(ind_peaks > i);
                tmpind2 = find(types_sorted(ind_edge_start(ind_peaks(tmpind))) == types_sorted(ind_edge_start(i)));
                if ~isempty(tmpind2)
                    peak_distance = [peak_distance xs(locs_sorted(ind_edge_start(ind_peaks(tmpind(tmpind2(1))))))-xs(locs_sorted(ind_edge_start(i)))];
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
        amp = [amp range_smoothed]; % Only amplitude
        
%         if isig == 80
%             close all
%             figure
%             subplot(1,3,[1 2])
%             plot(xs,c_smoothed_eval(:,isig))
%             hold on
%             plot(timestamp(range_ind),c_signal_woNharm(:,isig),'-x')
%             hold on
%             mycolor = lines(length(ispeak));
%             for i = 1:length(ind_edge_start)
%                 myrange = find(candidate_left(ind_edge_start(i))==c_smoothed_eval(:,isig)):find(candidate_right(ind_edge_end(i))==c_smoothed_eval(:,isig));
%                 if ispeak(i)
%                     plot(xs(myrange),c_smoothed_eval(myrange,isig),'Color',mycolor(i,:),'LineWidth',2)
%                 end
%             end
%             subplot(1,3,3)
%             set(gca,'Visible','off')
%             text(.1,.9,sprintf('nEdges = %g',nEdges(end)))
%             text(.1,.8,sprintf('SNR = %g',SNR(end)))
%             text(.1,.7,sprintf('Amplitude = %g',amp(end)))
%             text(.1,.6,sprintf('Peak dur mean = %g',peakdur_mean(end)))
%             text(.1,.5,sprintf('Peak dur std = %g',peakdur_std(end)))
%             text(.1,.4,sprintf('Peak dis mean = %g',peakdis_mean(end)))
%             text(.1,.3,sprintf('Peak dis std = %g',peakdis_std(end)))
%             1;
%         end

    end
    
    max_nEdges = nbasis - 5;
    max_SNR = 40;
    max_amp = 0.04;
    max_pw = 0.1;
    max_peakdur = 300;
    max_peakdis = 300;
    peakdur_mean(isnan(peakdur_mean)) = max_peakdur;
    peakdis_mean(isnan(peakdis_mean)) = max_peakdis;
    if ~log_trafo
        max_amp = 10.^max_amp;
    end
    nEdges(nEdges == 0) = .5;
    
    weight_edg = 1; % Weight between nEdges and SNR (Determine by PCA later?)
    weight_snr = 1;
    weight_amp = 1;
    weight_pw = 1;
    weight_peakdur = 1;
    weight_peakdis = 1;
    
    a = getcoef(smoothed_data_woNharm);
    b = a(2:end,:);
    
    pw = sqrt(sum((a(1:end-1,:)-b).^2,1));
    
%     princomp([log10(nEdges)' log10(SNR)' log10(amp)'])
    
%     radial_dist = (nEdges/(max_nEdges*1.05)).^weight_edg .* (SNR/(max_SNR*1.05)).^weight_snr .* (amp/(max_amp*1.05)).^weight_amp .* (pw/(max_pw*1.05)).^weight_pw .* (1./(peakdur_mean/(max_peakdur*1.05))).^weight_peakdur .* (1./(peakdis_mean/(max_peakdis*1.05))).^weight_peakdis;

    radial_dist = (nEdges/(max_nEdges*1.05)).^weight_edg .* (SNR/(max_SNR*1.05)).^weight_snr .* (amp/(max_amp*1.05)).^weight_amp .* (1./(peakdur_mean/(max_peakdur*1.05))).^weight_peakdur .* (1./(peakdis_mean/(max_peakdis*1.05))).^weight_peakdis;

    radial_dist = radial_dist.^(1/3); % To make distances smaller
    
end