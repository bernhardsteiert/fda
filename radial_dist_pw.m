function [radial_dist c_signal_woNharm range_ind] = radial_dist(isite,myextension,timeshift)
    if(~exist('myextension','var'))
        myextension = '';
    else
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
    
    a = getcoef(smoothed_data_woNharm);
    b = a(2:end,:);
    
    radial_dist = sqrt(sum((a(1:end-1,:)-b).^2,1));
%     radial_dist = sqrt(sum(getcoef(smoothed_data_woNharm).^2,1));

end