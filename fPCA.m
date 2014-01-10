function [scores c_signal range_ind] = fPCA(isite)
    
    load('harm_basis_fPCA.mat') % Contains only harm_basis_fPCA from all data-sets
    
    remotepath = mypath();
    
    warning('off','MATLAB:dispatcher:pathWarning')
    
    fdaMPath = [remotepath 'fda'];
    addpath(fdaMPath)

    grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
    addpath(grabdataPath)

    log_trafo = 1; % log-transform signal
    time_range = getbasisrange(harm_basis_fPCA);
    
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
    
    [tmp range_ind_min] = min(abs(timestamp - time_range(1)));
    [tmp range_ind_max] = min(abs(timestamp - time_range(2)));
    range_ind = range_ind_min:range_ind_max;
    
    if timestamp(range_ind(1)) < time_range(1)
        range_ind = range_ind(2:end);
    end
    if timestamp(range_ind(end)) > time_range(2)
        range_ind = range_ind(1:end-1);
    end
    
    smoothed_additional = smooth_basis(timestamp(range_ind),c_signal(range_ind,:),harm_basis_fPCA);
    
    scores = getcoef(smoothed_additional);

end