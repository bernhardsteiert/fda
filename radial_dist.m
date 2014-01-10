function [radial_dist c_signal_woNharm range_ind] = radial_dist(isite)
    
    load('harm_basis.mat') % Contains only harm_basis from all data-sets
    
    remotepath = mypath();
    
    warning('off','MATLAB:dispatcher:pathWarning')
    
    fdaMPath = [remotepath 'fda'];
    addpath(fdaMPath)

    grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
    addpath(grabdataPath)

    log_trafo = 1; % log-transform signal
    time_range = getbasisrange(harm_basis);
    
    if exist(remotepath,'dir')
        [timestamp,intensity] = grabdata(isite);
    else
        load(['./Workspaces/site_' num2str(isite)])
    end

    if log_trafo
        c_signal = log10(intensity);
    else
        c_signal = intensity;
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
    nbasis = 40;
    basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
    smoothed_data_woNharm = smooth_basis(timestamp(range_ind),c_signal_woNharm,basis);
    
    radial_dist = sqrt(sum(getcoef(smoothed_data_woNharm).^2,1));

end