function out_c_signal = register_signal(c_signal,extension)
% Step detection
% By using first derivative
% Criteria:
% - If 1 slope is 60% higher than all others --> Take this one
% - Define blocks by slope > -0.3 * max(slope) --> Small negative slopes are tolerated
% - If 1 block of rise is 2 timesteps or longer --> Take expectation of it (not sure about that ...)
% - If maximum is after minimum ... then do nothing? (not yet implemented)

if(~exist('extension','var'))
    out_c_signal = c_signal;
    return
elseif isempty(extension)
    out_c_signal = c_signal;
    return
end

switch extension
    case '12-25-2013'
        x = 45:55;
        % timestamp(50) = 128min = IC50 at harm_basis_fPCA
        ic50 = 50;
    case '12-08-2013'
        x = 2:12;
        ic50 = 6;
        % timestamp(6) = 128min = IC50 at harm_basis_fPCA
    otherwise
        out_c_signal = c_signal;
        return
end
        
out_c_signal = c_signal;

ind = 1:length(x)-1;
f = c_signal(x,:);
selected = zeros(size(f,2),1);
dx = x(2)-x(1);
newx = x(2:end)-dx/2;

dfdx = f(ind+1,:)-f(ind,:)./dx;
% plot(x,f)
% plot(newx,dfdx)

% return

fpos = dfdx;
maxcrit = abs(fpos) > .6*repmat(max(fpos,[],1),size(fpos,1),1);
selected(sum(maxcrit,1) == 1) = 1;
efposnorm = nan(size(selected));
tmpind = find(maxcrit(:,selected==1));
tmpx = repmat(newx',1,length(tmpind));
efposnorm(selected==1) = tmpx(tmpind);
% fpos(fpos<0) = 0; % All negative slopes are ignored; Too stringend!
fpos(fpos<-.3*repmat(max(fpos,[],1),size(fpos,1),1)) = nan; % Define blocks
% plot(newx,fpos)
% return
% fpos(fpos<.3*repmat(max(fpos,[],1),size(fpos,1),1)) = 0; % Second cut-off with threshold relative to maximum slope
fposnorm = fpos./repmat(nansum(fpos,1)*dx,size(fpos,1),1);

% fposnorm(fposnorm < .1) = 0; % Second cut-off after normalization with hard threshold

ind_to_test = find(selected==0)';
for i = ind_to_test
    
    blockstart = strfind([0 ~isnan(fposnorm(:,i)') 0], [0 1]);
    blockend = strfind([0 ~isnan(fposnorm(:,i)') 0], [1 0]);
    blocklength = blockend-blockstart;
    
    [maxlength blockpos] = max(blocklength);
    
    indstoremove = ones(size(fposnorm,1),1);
    if maxlength > 2
        selected(i) = 1;
        indstoremove(blockstart(blockpos):blockend(blockpos)-1) = 0;
    end
    fposnorm(indstoremove==1 | fposnorm(:,i)<0,i) = 0;
    
end

fposnorm = fposnorm./repmat(nansum(fposnorm,1)*dx,size(fpos,1),1);
% plot(newx,fposnorm)
% return
efposnorm(ind_to_test) = nansum(fposnorm(:,ind_to_test) .* repmat(newx',1,size(fposnorm(:,ind_to_test),2)) * dx,1);

% hold on
% ylim = get(gca,'YLim');
% 
% colors = lines(length(efposnorm));
% for i = 1:length(efposnorm)
%     
%     plot([efposnorm(i) efposnorm(i)],ylim,'--','Color',colors(i,:))
%     
% end

% close all
nsubplots = size(f,2);
% nrows = ceil(nsubplots^rowstocols);
% ncols = ceil(nsubplots / nrows);

for i = 1:nsubplots
%     subplot(nrows,ncols,i)
%     plot(x,f(:,i))
%     hold on
%     ylim = get(gca,'YLim');

%     if max((strfind([0 (fposnorm(:,i)' > 0) 0], [1 0]) - strfind([0 (fposnorm(:,i)' > 0) 0], [0 1]))) > 1
    if selected(i)
%         plot([efposnorm(i) efposnorm(i)],ylim,'--')
        deltaind = -(ic50 - round(efposnorm(i)));
        if deltaind > 0
            out_c_signal(1:size(c_signal,1)-deltaind,i) = c_signal(deltaind+1:size(c_signal,1),i);
        elseif deltaind < 0
            out_c_signal(-deltaind+1:size(c_signal,1),i) = c_signal(1:size(c_signal,1)+deltaind,i);
    end
end

end