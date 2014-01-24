function [radial_dist c_signal_woNharm range_ind] = freq_analysis(isite,freq,myextension,timeshift)
    if(~exist('myextension','var'))
        freq = 2e-4;
    end
    if(~exist('myextension','var'))
        myextension = '';
    else
        myextension = ['_' myextension];
    end
    if(~exist('timeshift','var'))
        timeshift = 0;
    end
    
    remotepath = mypath();
    
    warning('off','MATLAB:dispatcher:pathWarning')
    
    fdaMPath = [remotepath 'fda'];
    addpath(fdaMPath)

    grabdataPath = [remotepath 'Code + Stage and Outputsignal'];
    addpath(grabdataPath)

    log_trafo = 1; % log-transform signal
    register = 1; % register IC50
    time_range = [200 510];
    
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
    c_signal = c_signal - repmat(nanmean(c_signal,2),1,size(c_signal,2)); % Mean subtraction in y-direction
    c_signal_woNharm = c_signal - repmat(nanmean(c_signal,1),size(c_signal,1),1); % Mean subtraction in x-direction
    
    [tmp range_ind_min] = min(abs(timestamp - time_range(1)));
    [tmp range_ind_max] = min(abs(timestamp - time_range(2)));
    range_ind = range_ind_min:range_ind_max;
    
    if timestamp(range_ind(1)) < time_range(1)
        range_ind = range_ind(2:end);
    end
    if timestamp(range_ind(end)) > time_range(2)
        range_ind = range_ind(1:end-1);
    end
    
    Fs = 1./((timestamp(2)-timestamp(1))*60); % Sampling every 5 min
    L = length(range_ind);

    NFFT = 2^nextpow2(L);
        
    Y = fft(c_signal_woNharm(range_ind,:),NFFT,1)/L;

    f = Fs/2*linspace(0,1,NFFT/2+1);
    [tmp indf] = min(abs(f-freq));
    Yabs = abs(Y);
    Yabs(Yabs < 5*nanmean(nanmean(Yabs(f > 1e-3,:)))) = 0;
    Yabs([find(f < 1.1e-4) find(f > 1e-3)],:) = 10*Yabs([find(f < 1.1e-4) find(f > 1e-3)],:); % Penalize slow & fast components
    radial_dist = 2*Yabs(indf,:)./sum(Yabs(setdiff(1:size(Yabs,1),indf),:),1);

end