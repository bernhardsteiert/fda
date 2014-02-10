function [scores c_signal range_ind] = fPCA_late(isite,myextension,timeshift)
    if(~exist('myextension','var'))
        myextension = '';
    elseif(~isempty(myextension))
        myextension = ['_' myextension];
    end
    if(~exist('timeshift','var'))
        timeshift = 0;
    end

    load('harm_basis.mat') % Contains only harm_basis_fPCA from all data-sets
    
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
    
    scores = getcoef(smoothed_additional);

end