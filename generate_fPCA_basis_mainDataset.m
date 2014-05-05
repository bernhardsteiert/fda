close all
clear all
clc

nharm = 5;

sites = [4	5	6	7	8	9,...
17	16	15	14	13	12,...
24	25	26	27	28	29,...
37	36	35	34	33	32,...
44	45	46	47	48	49,...
57	56	55	54	53	52,...
64	65	66	67	68	69]; % 
timeshift = 0;
flipharm = ones(1,nharm);
Rmat = diag(ones(1,nharm));

load(sprintf('FOXO3a_figures/Workspaces/site_%i.mat',sites(1)))
timestamp = timestamp - timeshift; % Shift to main data set

remotepath = mypath();

fdaMPath = [remotepath 'fda'];
addpath(fdaMPath)


grabdataPath = 'C:\Users\SS240\Dropbox (Somponnat workspace)\Bernhard-Pat H5 file sharing\130722';


times = cell(0);
signals = cell(0);
celltype = [];

for isite = sites
    if exist(grabdataPath,'dir')
        [timestamp,intensity] = grabdata_new(isite);
    else
        load(['./Workspaces/site_' num2str(isite)])
        times{end+1} = timestamp;
    end

    log_trafo = 1; % log-transform signal

    if log_trafo
        signals{end+1} = log10(intensity);
    else
        signals{end+1} = intensity;
    end
    
    celltype = [celltype ones(1,size(intensity,2))*isite];
end

c_signal = cell2mat(signals);
c_signal(isinf(c_signal)) = nan;

for i = 1:size(c_signal,2)
    c_signal(:,i) = interp1(timestamp(~isnan(c_signal(:,i))),c_signal(~isnan(c_signal(:,i)),i),timestamp);
    vec = ~isnan(c_signal(:,i))';
    rl = find(vec ~= [vec(2:end), vec(end)+1]);
    data =  vec(rl);
    rl(2:end) = rl(2:end) - rl(1:end-1);
    if ~data(1)
        c_signal(1:rl(1),i) = c_signal(rl(1)+1,i);
    end
    if ~data(end)
        c_signal(end-rl(end)+1:end,i) = c_signal(end-rl(end),i);
    end
end

time_range = [50 220];


[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

nbasis = round(length(range_ind)/1.5);

basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
smoothed_data = smooth_basis(timestamp(range_ind),c_signal(range_ind,:),basis);

c_signal_pcastr = pca_fd(smoothed_data, nharm, fdPar(basis, int2Lfd(2), 0), 0); % WITHOUT CENTERING!!

fdobj = getcoef(c_signal_pcastr.harmfd);
fdobj = repmat(flipharm,size(fdobj,1),1).*fdobj;
basisobj = getbasis(c_signal_pcastr.harmfd);

newcoef = fdobj;
for irow = 1:size(Rmat,1)
    newcoef(:,irow) = sum(repmat(Rmat(irow,:),size(fdobj,1),1) .* fdobj(:,1:size(Rmat,1)),2);
end

harm_fPCA = fd(newcoef,basisobj,{'time','Rotated harmonics','Rotated harmonics for variables'});

harm_basis = create_fd_basis(harm_fPCA);
plot(harm_basis)

myextension = '';
% return
save(sprintf('harm_basis_%s',myextension),'harm_basis')