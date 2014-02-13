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
    case '01-22-2014'
        x = 7:17;
        ic50 = 12;
        % timestamp(12) = 128min = IC50 at harm_basis_fPCA
    case '01-26-2014'
        x = 59:69;
        ic50 = 64;
        % timestamp(65) = 128min = IC50 at harm_basis_fPCA
    case '01-27-2014'
        x = 91:101;
        ic50 = 96;
        % timestamp(96) = 128min = IC50 at harm_basis_fPCA
    otherwise
        out_c_signal = c_signal;
        return
end
        
out_c_signal = c_signal;

ind = 1:length(x)-1;
f = c_signal(x,:);
selected = zeros(size(f,2),1);
slope = nan(size(f,2),1);
dx = x(2)-x(1);
newx = x(2:end)-dx/2;

dfdx = f(ind+1,:)-f(ind,:)./dx;
% plot(x,f)
% plot(newx,dfdx)

% return

fpos = dfdx;
% maxcrit = abs(fpos) > .6*repmat(max(fpos,[],1),size(fpos,1),1);
maxcrit = fpos > .6*repmat(max(fpos,[],1),size(fpos,1),1);
selected(sum(maxcrit,1) == 1) = 1;
slope(selected==1) = max(fpos(:,selected==1),[],1);
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
    if maxlength > 2
        blockstart = strfind([0 fposnorm(:,i)'>0 0], [0 1]);
        blockend = strfind([0 fposnorm(:,i)'>0 0], [1 0]);

        if ~isempty(blockstart)
            slope(i) = (f(blockend(end)-1) - f(blockstart(1)))./(dx*(blockend(end)-blockstart(1)));
            if slope(i) < 0
                slope(i) = nan;
                selected(i) = 0;
            end
        else
            selected(i) = 0;
        end
    end
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
% figure
nsubplots = size(f,2);
rowstocols = .5;
nrows = ceil(nsubplots^rowstocols);
ncols = ceil(nsubplots / nrows);
deltaind = zeros(size(selected));

for i = 1:nsubplots
%     subplot(nrows,ncols,i)
%     plot(x,f(:,i))
%     hold on
%     ylim = get(gca,'YLim');

%     if max((strfind([0 (fposnorm(:,i)' > 0) 0], [1 0]) - strfind([0 (fposnorm(:,i)' > 0) 0], [0 1]))) > 1
    if selected(i)
%         plot([efposnorm(i) efposnorm(i)],ylim,'--')
%         plot([round(efposnorm(i)-1e-10) round(efposnorm(i)-1e-10)],ylim)
%         plot([ic50 ic50],ylim,'r--')
        deltaind(i) = -(ic50 - round(efposnorm(i)-1e-10));
        if deltaind(i) > 0
            out_c_signal(1:size(c_signal,1)-deltaind(i),i) = c_signal(deltaind(i)+1:size(c_signal,1),i);
            % Negative time-shift data looks always ugly --> artifact
%         elseif deltaind(i) < 0
%             out_c_signal(-deltaind(i)+1:size(c_signal,1),i) = c_signal(1:size(c_signal,1)+deltaind(i),i);
        end
%         plot(x,c_signal(x+deltaind(i),i),'r')
    end

end

return

figure
unidelta = setdiff(unique(deltaind),0);
slopemat = [];
for u = 1:length(unidelta)
    slopemat = padconcatenation(slopemat,slope(deltaind == unidelta(u)),2);
end

plot(deltaind,slope,'*')
% plot(deltaind(deltaind~=0),slope(deltaind~=0),'*')

% sum(selected)
% find(deltaind~=0)
% boxplot(slopemat)
% set(gca,'XTick',1:length(unidelta),'XTickLabel',unidelta)