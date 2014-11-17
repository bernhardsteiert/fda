%% Dual-sensor first steps
% load('c_signal_dual')
% save('c_signal_dual_2','celltype','timestamp','c_signal_foxo','c_signal_ekarev')
load('c_signal_dual_2')

exten = '02-15-2014';
obs = {'Dual_FOXO3a','Dual_EKAREV'};
sites = [7	8	9   12	11	10  25	26	27  30	29	28  43	44	45  48	47	46  61	62	63];

timeshift = -72;
timestamp = timestamp - timeshift; % Shift to main data set

% Remove Infs / Zeros
c_signal_foxo(isinf(c_signal_foxo)) = nan;
c_signal_ekarev(isinf(c_signal_ekarev)) = nan;
% c_signal_ekarev(c_signal_ekarev==0) = nan;

% Normalization
time_cutoff = 200;
myind = timestamp>time_cutoff;
c_signal_foxo = c_signal_foxo./repmat(range(c_signal_foxo(myind,:),1),size(c_signal_foxo,1),1);
c_signal_foxo = c_signal_foxo - repmat(min(c_signal_foxo(myind,:),[],1),size(c_signal_foxo,1),1);
c_signal_ekarev = c_signal_ekarev./repmat(range(c_signal_ekarev(myind,:),1),size(c_signal_ekarev,1),1);
c_signal_ekarev = c_signal_ekarev - repmat(min(c_signal_ekarev(myind,:),[],1),size(c_signal_ekarev,1),1);

%% Plot some traces (BTC)
f1 = figure;

xfac = 1.5;
yfac = 1;
fontsize = 10;

setFigure(f1,xfac,yfac,fontsize)

nrow = 4;
ncol = 5;

myind2 = 1:length(celltype);
myind2 = find(celltype == 9); % Only plot BTC 100

for irow = 1:nrow
    for icol = 1:ncol
        subplot(nrow,ncol,(irow-1)*ncol+icol)
        hold on
        
        plot(timestamp,c_signal_foxo(:,myind2((irow-1)*ncol+icol)),'b')
        plot(timestamp,c_signal_ekarev(:,myind2((irow-1)*ncol+icol)),'r')
%         plot(c_signal_foxo(:,myind((irow-1)*ncol+icol)),c_signal_ekarev(:,myind((irow-1)*ncol+icol)))

        set(gca,'XLim',[200 1500])
    end
end

%% Outlier removal

c_signal_foxo(isinf(c_signal_foxo)) = nan;
removed_inds = nan(size(c_signal_foxo));

c_signal_foxo_interp = c_signal_foxo;
for i = 1:size(c_signal_foxo,2)
    if sum(~isnan(c_signal_foxo_interp(:,i))) > 1
        c_signal_foxo_interp(:,i) = interp1(timestamp(~isnan(c_signal_foxo(:,i))),c_signal_foxo(~isnan(c_signal_foxo(:,i)),i),timestamp);
        vec = ~isnan(c_signal_foxo_interp(:,i))';
        rl = find(vec ~= [vec(2:end), vec(end)+1]);
        data =  vec(rl);
        rl(2:end) = rl(2:end) - rl(1:end-1);
        if ~data(1)
            c_signal_foxo_interp(1:rl(1),i) = c_signal_foxo_interp(rl(1)+1,i);
        end
        if ~data(end)
            c_signal_foxo_interp(end-rl(end)+1:end,i) = c_signal_foxo_interp(end-rl(end),i);
        end

        % Outlier correction
        mysignal = c_signal_foxo_interp(:,i);
        testdiff = diff(c_signal_foxo_interp(:,i));

        win = 15;

        allranges = nan(1,length(testdiff)-win);
        for islide = 1:length(testdiff)-win
            allranges(islide) = range(testdiff(islide:islide+win-1));
        end

        range_sorted = sort(allranges);
        range_thres = range_sorted(round(length(range_sorted)*.5))-min(range_sorted);

        outliers = find(allranges > min(range_sorted)+4*range_thres);
        cleaned_inds = ones(1,length(mysignal));

        if ~isempty(outliers)
            outlier_bd = find(diff(outliers) > 1);

            outlier_lb = [outliers(1) outliers(outlier_bd+1)];
            outlier_ub = [outliers(outlier_bd) outliers(end)];

            for iout = 1:length(outlier_lb)
                if outlier_lb(iout) > win/2
                    cleaned_inds(round((outlier_lb(iout):outlier_ub(iout))+win/2)) = 0;
                else
                    cleaned_inds(outlier_lb(iout):round(outlier_ub(iout)+win/2)) = 0;
                end
            end
        end
        cleaned_inds = logical(cleaned_inds);
        cleaned_inds(timestamp < 200) = 1;
        c_signal_foxo(~cleaned_inds,i) = nan;
        removed_inds(:,i) = ~cleaned_inds;
    end
end

onlyNaN = sum(~isnan(c_signal_foxo),1) == 0 | sum(~isnan(c_signal_ekarev),1) == 0;
c_signal_foxo = c_signal_foxo(:,~onlyNaN);
c_signal_ekarev = c_signal_ekarev(:,~onlyNaN);
celltype = celltype(~onlyNaN);
removed_inds = removed_inds(:,~onlyNaN);

%% Interpolate between missing points
% Neccessary for spline interpolation to work
% For later analysis, interpolated data can be neglected again
myind = timestamp>time_cutoff;

c_signal_ekarev_interp = c_signal_ekarev;
c_signal_foxo_interp = c_signal_foxo;
removedInds = [];
for i = 1:size(c_signal_ekarev,2)
    remove = 0;
    if sum(~isnan(c_signal_ekarev(myind,i))) > 2 && range(c_signal_ekarev(myind,i)) ~= 0
        c_signal_ekarev_interp(:,i) = interp1(timestamp(~isnan(c_signal_ekarev(:,i))),c_signal_ekarev(~isnan(c_signal_ekarev(:,i)),i),timestamp);
        vec = ~isnan(c_signal_ekarev_interp(:,i))';
        rl = find(vec ~= [vec(2:end), vec(end)+1]);
        data =  vec(rl);
        rl(2:end) = rl(2:end) - rl(1:end-1);
        if ~data(1)
            c_signal_ekarev_interp(1:rl(1),i) = c_signal_ekarev_interp(rl(1)+1,i);
        end
        if ~data(end)
            c_signal_ekarev_interp(end-rl(end)+1:end,i) = c_signal_ekarev_interp(end-rl(end),i);
        end
    else
        remove = 1;
    end
    if sum(~isnan(c_signal_foxo(myind,i))) > 2 && range(c_signal_foxo(myind,i)) ~= 0
        c_signal_foxo_interp(:,i) = interp1(timestamp(~isnan(c_signal_foxo(:,i))),c_signal_foxo(~isnan(c_signal_foxo(:,i)),i),timestamp);
        vec = ~isnan(c_signal_foxo_interp(:,i))';
        rl = find(vec ~= [vec(2:end), vec(end)+1]);
        data =  vec(rl);
        rl(2:end) = rl(2:end) - rl(1:end-1);
        if ~data(1)
            c_signal_foxo_interp(1:rl(1),i) = c_signal_foxo_interp(rl(1)+1,i);
        end
        if ~data(end)
            c_signal_foxo_interp(end-rl(end)+1:end,i) = c_signal_foxo_interp(end-rl(end),i);
        end
    else
        remove = 1;
    end
    
    if remove
        removedInds = [removedInds i];
    end
end

keepInds = setdiff(1:size(c_signal_foxo,2),removedInds);
c_signal_foxo = c_signal_foxo(:,keepInds);
c_signal_foxo_interp = c_signal_foxo_interp(:,keepInds);
c_signal_ekarev = c_signal_ekarev(:,keepInds);
c_signal_ekarev_interp = c_signal_ekarev_interp(:,keepInds);
celltype = celltype(keepInds);
removed_inds = removed_inds(:,keepInds);

c_signal_foxo_interp = c_signal_foxo_interp./repmat(range(c_signal_foxo_interp(myind,:),1),size(c_signal_foxo_interp,1),1);
c_signal_foxo_interp = c_signal_foxo_interp - repmat(min(c_signal_foxo_interp(myind,:),[],1),size(c_signal_foxo_interp,1),1);
c_signal_ekarev_interp = c_signal_ekarev_interp./repmat(range(c_signal_ekarev_interp(myind,:),1),size(c_signal_ekarev_interp,1),1);
c_signal_ekarev_interp = c_signal_ekarev_interp - repmat(min(c_signal_ekarev_interp(myind,:),[],1),size(c_signal_ekarev_interp,1),1);

%% Plot some interpolated traces
f2 = figure;

xfac = 1.5;
yfac = 1;
fontsize = 10;

setFigure(f2,xfac,yfac,fontsize)

nrow = 4;
ncol = 5;

myind = 1:length(celltype);
myind = find(celltype == 9); % Only plot BTC 100

for irow = 1:nrow
    for icol = 1:ncol
        subplot(nrow,ncol,(irow-1)*ncol+icol)
        hold on
        
        plot(timestamp,c_signal_foxo_interp(:,myind((irow-1)*ncol+icol)),'b')
        plot(timestamp,c_signal_ekarev_interp(:,myind((irow-1)*ncol+icol)),'r')
%         plot(c_signal_foxo(:,myind((irow-1)*ncol+icol)),c_signal_ekarev(:,myind((irow-1)*ncol+icol)))

        set(gca,'XLim',[200 1500])
    end
end

%% fPCA basis calculation
flipharm = [1 1 1 1 1];
% Rmat = diag(flipharm);
time_range = [200 timestamp(end)];
[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

nbasis = round(length(range_ind)/1.5);

basis = create_bspline_basis([timestamp(range_ind(1)) timestamp(range_ind(end))], nbasis);
smoothed_data_foxo = smooth_basis(timestamp(range_ind),c_signal_foxo_interp(range_ind,:),basis);
smoothed_data_ekarev = smooth_basis(timestamp(range_ind),c_signal_ekarev_interp(range_ind,:),basis);

nharm = length(flipharm);
c_signal_foxo_pcastr = pca_fd(smoothed_data_foxo, nharm, fdPar(basis, int2Lfd(2), 0), 0); % WITHOUT CENTERING!!
c_signal_ekarev_pcastr = pca_fd(smoothed_data_ekarev, nharm, fdPar(basis, int2Lfd(2), 0), 0); % WITHOUT CENTERING!!

fdobj_foxo = getcoef(c_signal_foxo_pcastr.harmfd);
fdobj_foxo = repmat(flipharm,size(fdobj_foxo,1),1).*fdobj_foxo;
basisobj_foxo = getbasis(c_signal_foxo_pcastr.harmfd);
fdobj_ekarev = getcoef(c_signal_ekarev_pcastr.harmfd);
fdobj_ekarev = repmat(flipharm,size(fdobj_ekarev,1),1).*fdobj_ekarev;
basisobj_ekarev = getbasis(c_signal_ekarev_pcastr.harmfd);

newcoef_foxo = fdobj_foxo;
newcoef_ekarev = fdobj_ekarev;
% for irow = 1:size(Rmat,1)
%     newcoef(:,irow) = sum(repmat(Rmat(irow,:),size(fdobj,1),1) .* fdobj(:,1:size(Rmat,1)),2);
% end

harm_fPCA_foxo = fd(newcoef_foxo,basisobj_foxo,{'time','Harmonics','Harmonics for variables'});
harm_fPCA_ekarev = fd(newcoef_ekarev,basisobj_ekarev,{'time','Harmonics','Harmonics for variables'});

harm_basis_foxo = create_fd_basis(harm_fPCA_foxo);
harm_basis_ekarev = create_fd_basis(harm_fPCA_ekarev);
figure
plot(harm_basis_foxo)
figure
plot(harm_basis_ekarev)

%% fPCA basis removal
smoothed_additional_foxo = smooth_basis(timestamp(range_ind),c_signal_foxo_interp(range_ind,:),harm_basis_foxo);

harm_eval_foxo = eval_basis(harm_basis_foxo,timestamp(range_ind));

fitcoef_foxo = getcoef(smoothed_additional_foxo);
data_fpca_repr_foxo = fitcoef_foxo'*harm_eval_foxo';

c_signal_foxo_woNharm = c_signal_foxo_interp(range_ind,:)-data_fpca_repr_foxo';
c_signal_foxo_woNharm_NaN = c_signal_foxo_woNharm;
c_signal_foxo_woNharm_NaN(logical(removed_inds(range_ind,:))) = nan;


smoothed_additional_ekarev = smooth_basis(timestamp(range_ind),c_signal_ekarev_interp(range_ind,:),harm_basis_ekarev);

harm_eval_ekarev = eval_basis(harm_basis_ekarev,timestamp(range_ind));

fitcoef_ekarev = getcoef(smoothed_additional_ekarev);
data_fpca_repr_ekarev = fitcoef_ekarev'*harm_eval_ekarev';

c_signal_ekarev_woNharm = c_signal_ekarev_interp(range_ind,:)-data_fpca_repr_ekarev';
c_signal_ekarev_woNharm_NaN = c_signal_ekarev_woNharm;
c_signal_ekarev_woNharm_NaN(logical(removed_inds(range_ind,:))) = nan;

%% Plot some fPCA basis cleaned traces (BTC)
f2 = figure;

xfac = 1.5;
yfac = 1;
fontsize = 10;

setFigure(f2,xfac,yfac,fontsize)

nrow = 4;
ncol = 5;

myind = 1:length(celltype);
myind = find(celltype == 9); % Only plot BTC 100

for irow = 1:nrow
    for icol = 1:ncol
        subplot(nrow,ncol,(irow-1)*ncol+icol)
        hold on
        
        plot(timestamp(range_ind),c_signal_foxo_woNharm_NaN(:,myind((irow-1)*ncol+icol)),'b')
        plot(timestamp(range_ind),c_signal_ekarev_woNharm_NaN(:,myind((irow-1)*ncol+icol)),'r')
        set(gca,'XLim',[200 1500])
    end
end


f2 = figure;

setFigure(f2,xfac,yfac,fontsize)

for irow = 1:nrow
    for icol = 1:ncol
        subplot(nrow,ncol,(irow-1)*ncol+icol)
        hold on
        
        plot(c_signal_foxo_woNharm_NaN(:,myind((irow-1)*ncol+icol)),c_signal_ekarev_woNharm_NaN(:,myind((irow-1)*ncol+icol)),'k*')
        mycor = corrcoef(c_signal_foxo_woNharm_NaN(~isnan(c_signal_foxo_woNharm_NaN(:,myind((irow-1)*ncol+icol))),myind((irow-1)*ncol+icol)),c_signal_ekarev_woNharm_NaN(~isnan(c_signal_ekarev_woNharm_NaN(:,myind((irow-1)*ncol+icol))),myind((irow-1)*ncol+icol)));
        text(.5,.8,sprintf('CorrCoef = %g',mycor(1,2)),'Units','normalized')

    end
end

%% Analyze correlation between EKAREV and FOXO3a
myCorrC = nan(1,size(c_signal_foxo_woNharm_NaN,2));
for i = 1:length(myCorrC)
    mycor = corrcoef(c_signal_foxo_woNharm_NaN(~isnan(c_signal_foxo_woNharm_NaN(:,i)),i),c_signal_ekarev_woNharm_NaN(~isnan(c_signal_ekarev_woNharm_NaN(:,i)),i));
    myCorrC(i) = mycor(1,2);
end
[tmp ind_corrC] = sort(myCorrC,'descend');

%% Plot some fPCA basis cleaned traces (most correlated)
f2 = figure;

xfac = 1.5;
yfac = 1;
fontsize = 10;

setFigure(f2,xfac,yfac,fontsize)

nrow = 4;
ncol = 5;

myind = ind_corrC; % Most correlated traces

for irow = 1:nrow
    for icol = 1:ncol
        subplot(nrow,ncol,(irow-1)*ncol+icol)
        hold on
        
        plot(timestamp(range_ind),c_signal_foxo_woNharm_NaN(:,myind((irow-1)*ncol+icol)),'b')
        plot(timestamp(range_ind),c_signal_ekarev_woNharm_NaN(:,myind((irow-1)*ncol+icol)),'r')
        set(gca,'XLim',[200 1500])
        
%         plot(c_signal_foxo_woNharm_NaN(:,myind((irow-1)*ncol+icol)),c_signal_ekarev_woNharm_NaN(:,myind((irow-1)*ncol+icol)),'k*')
%         mycor = corrcoef(c_signal_foxo_woNharm_NaN(:,myind((irow-1)*ncol+icol)),c_signal_ekarev_woNharm_NaN(:,myind((irow-1)*ncol+icol)));
%         text(.5,.8,sprintf('CorrCoef = %g',mycor(1,2)),'Units','normalized')
    end
end

%% Boxplot for high doses
highdoses = [9 10 27 28 45 46 63];
myext = '02-15-2014_retracked';

mylab = {};
boxCorrC = [];

for i = 1:length(highdoses)
    s = siteprop(highdoses(i),myext);
    
    mylab{end+1} = s.lig_name;
    boxCorrC = padconcatenation(boxCorrC,myCorrC(celltype == highdoses(i)),1);
    
end

figure
boxplot(boxCorrC')

set(gca,'XTick',1:length(highdoses),'XTickLabel',mylab)