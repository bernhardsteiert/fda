%% Immunostaining data of fixed parental cells
% cell_name = {'184A1','BT20','HCC1806','HS578T','MCF7','MCF10A','MDA231','SKBR3','T47D'}; % FOXO3a of MDA231 seems not to be normalized
cell_name = {'184A1','BT20','HCC1806','HS578T','MCF7','MCF10A12202013','SKBR3','T47D'};
% cell_name = {'184A1','BT20','HCC1806','HS578T','MCF7','MCF10A','MCF10A12202013','SKBR3','T47D'};

Ligand = {'EGF';'IGF1';'FGF1';'HRG';'HGF';'EPR';'BTC';'NS'}; %row 1-8
Tx = {'No Drug';'AKTi';'MEKi';'Both'}; % col 1-4
obs = {'FOXO3a', 'pERK', 'pAKT'};
% obs = {'pERK', 'pAKT'};
timepoints_184a1 = [0, 5, 10, 15, 20, 30, 45, 60, 90, 120, 180, 300, 480];
timepoints_all   = [0,        15,     30,     60, 90, 120, 180, 240];

close all

nColor = length(cell_name);

nrows = length(Ligand);
ncols = length(Tx);

levels = cell(length(cell_name),length(obs));
alliqrs = [];

for ic = 1:length(cell_name)
    
    for iobs = 1:length(obs)

    levels{ic,iobs} = nan(nrows,ncols,3); % 3 levels: t(0); t(15); mean(t(90:end))

        if ic == 1
            timepoints = timepoints_184a1;
            tlt80 = 8;
        else
            timepoints = timepoints_all;
            tlt80 = 4;
        end

        medians = nan(length(timepoints),nrows,ncols);
        iqrs = nan(length(timepoints),nrows,ncols);

        data = load(['./Workspaces/' cell_name{ic} '_' obs{iobs}]);
        mydata = getfield(data,['single_' obs{iobs}]);

        for ilig = 1:nrows
            for idrug = 1:ncols
                for it = 1:length(timepoints)
                    data = mydata{it,ilig,idrug}(:,1);
                    medians(it,ilig,idrug) = median(10.^data);
                    iqrs(it,ilig,idrug) = iqr(10.^data);
                end
            end
        end
        
        if strmatch(obs{iobs},'pERK')
            % pERK
            medians = medians - mean(medians(1,:,3)); % timepoint 1; all Ligands; MEKi
            medians = medians ./ max(medians(:));
        elseif strmatch(obs{iobs},'pAKT')
            % pAKT
            medians = medians - mean(medians(1,:,2)); % timepoint 1; all Ligands; MEKi
            medians = medians ./ max(medians(:));
        elseif strmatch(obs{iobs},'FOXO3a')
            % FOXO3a
            rangefoxo = range(medians(:));
            iqrs = iqrs ./ rangefoxo;
            medians = medians ./ rangefoxo;
        end

        for ilig = 1:nrows
            for idrug = 1:ncols
                
                levels{ic,iobs}(ilig,idrug,1) = medians(1,ilig,idrug);
                levels{ic,iobs}(ilig,idrug,2) = medians(timepoints == 15,ilig,idrug);
                levels{ic,iobs}(ilig,idrug,3) = mean(medians(timepoints > 80,ilig,idrug));
                if strmatch(obs{iobs},'FOXO3a') 
                    alliqrs = [alliqrs mean(iqrs(timepoints > 80,ilig,idrug))];
                end
                
            end
        end
        
    end
end


alllevels = nan(length(cell_name)*length(Ligand)*length(Tx),3,length(obs));
icount = 1;
for iobs = 1:length(obs);
    cell_inds = [];
    lig_inds = [];
    drug_inds = [];
    for il = 1:size(levels,1)
        mylevel = levels{il,iobs};
        for i = 1:size(mylevel,1)
            for j = 1:size(mylevel,2)
                tmplevel = [];
                for k = 1:size(mylevel,3)
                   alllevels(icount,k,iobs) = mylevel(i,j,k);
                end
                icount = icount + 1;
                
                cell_inds = [cell_inds il];
                lig_inds = [lig_inds i];
                drug_inds = [drug_inds j];
            end
        end
    end
    icount = 1;
end

% ------------------------------------

% FOXO3a from ERK & AKT

% Choose FOXO3a IQRs that are regressed
regressiqrs = [1 3]; % 1: WT; 2: AKTi; 3: MEKi; 4: Both

iqrinds = 1:4:length(alliqrs);
choseniqrs = [];
for iregiqr = regressiqrs
    choseniqrs = [choseniqrs iqrinds+iregiqr-1];
end

% Choose data that is used for regression

% Use both AKT + ERK data for regression
% regressvecs = [ones(size(alllevels,1),1) alllevels(:,:,2) alllevels(:,:,3)];
% regressvecs_iqr = [ones(length([iqrinds iqrinds+2]),1) alllevels([iqrinds iqrinds+2],:,2) alllevels([iqrinds iqrinds+2],:,3)];
% titstr = 'ERK + AKT';

% Use only AKT data for regression
regressvecs = [ones(size(alllevels,1),1) alllevels(:,:,3)];
regressvecs_iqr = [ones(length(regressiqrs)*length(iqrinds),1) alllevels(choseniqrs,:,3)];
titstr = 'Only AKT';

% Use only ERK data for regression
% regressvecs = [ones(size(alllevels,1),1) alllevels(:,:,2)];
% regressvecs_iqr = [ones(length([iqrinds iqrinds+2]),1) alllevels([iqrinds iqrinds+2],:,2)];
% titstr = 'Only ERK';

par = [];
for ilevel = 1:size(alllevels,2)
    par = [par regress(alllevels(:,ilevel,1),regressvecs)];
end

par = [par regress(alliqrs(choseniqrs)',regressvecs_iqr)]; % only MEKi

ylabels = {'FOXO3a baseline','FOXO3a peakheight','FOXO3a steady-state','FOXO3a late IQR'};
for iplot = 1:size(par,2)
    subplot(2,2,iplot)
    hold on

    tmpx = regressvecs*par(:,iplot);
    if iplot < 4
        plot(tmpx,alllevels(:,iplot,1),'kx')
    end
    plot(1.1*[0 max(tmpx)],1.1*[0 max(tmpx)],'k:')
    
    axis equal
    title(titstr)
    ylabel(ylabels{iplot})
    xlabel('predictor')
end

% Color by ligand:
% colmap = hsv(length(Ligand));
% legh = [];
% for ilig = 1:length(Ligand)
%     plot(tmpx(lig_inds == ilig),alliqrs(lig_inds == ilig),'k.')
%     legh = [legh plot(tmpx(lig_inds == ilig & drug_inds == 1),alliqrs(lig_inds == ilig & drug_inds == 1),'o','MarkerFaceColor',colmap(ilig,:))]; % WT
%     plot(tmpx(lig_inds == ilig & drug_inds == 3),alliqrs(lig_inds == ilig & drug_inds == 3),'s','MarkerFaceColor',colmap(ilig,:))
% %     plot(tmpx(lig_inds == ilig & drug_inds == 2),alliqrs(lig_inds == ilig & drug_inds == 2),'d','MarkerFaceColor',colmap(ilig,:))
% end
% legend(legh,Ligand)

% Color by celltype:
colmap = hsv(length(cell_name));
legh = [];
for ic = 1:length(cell_name)
    plot(tmpx(cell_inds == ic),alliqrs(cell_inds == ic),'k.')
    legh = [legh plot(tmpx(cell_inds == ic & drug_inds == 1),alliqrs(cell_inds == ic & drug_inds == 1),'o','MarkerFaceColor',colmap(ic,:))]; % WT
    plot(tmpx(cell_inds == ic & drug_inds == 3),alliqrs(cell_inds == ic & drug_inds == 3),'s','MarkerFaceColor',colmap(ic,:))
%     plot(tmpx(cell_inds == ic & drug_inds == 2),alliqrs(cell_inds == ic & drug_inds == 2),'d','MarkerFaceColor',colmap(ic,:))
end
legend(legh,cell_name)


% Scoring function: (from diagonal)
% 0deg optimal; 90deg neutral; 180deg worst case
% score: sum(angle_i*length_i) --> only angles that have length are scored

ligscore = [];
for ilig = 1:length(Ligand)
    
    % Calculate and score angle and length between:
    p1x = tmpx(lig_inds == ilig & drug_inds == 1); % x of point 1
    p1y = alliqrs(lig_inds == ilig & drug_inds == 1)'; % y of point 1
    
    p2x = tmpx(lig_inds == ilig & drug_inds == 3); % x of point 2
    p2y = alliqrs(lig_inds == ilig & drug_inds == 3)'; % y of point 2
    
    optdir = 1./sqrt(2).*[-1; -1];
    ligdir = [p2x-p1x p2y-p1y];
    score_lengths = diag(ligdir*repmat(optdir,1,size(ligdir,1))); % positive for angles < 90; negative otherwise
    score_angles = pi-acos(score_lengths./(sqrt(ligdir(:,1).^2+ligdir(:,2).^2)));
    
    ligscore = [ligscore sum(score_lengths.*score_angles)];
    
end

figure
bar(ligscore)
set(gca,'XTick',1:length(Ligand),'XTickLabel',Ligand)
ylabel('Diagonal score')

ctscore = [];
for ic = 1:length(cell_name)
    
    % Calculate and score angle and length between:
    p1x = tmpx(cell_inds == ic & drug_inds == 1); % x of point 1
    p1y = alliqrs(cell_inds == ic & drug_inds == 1)'; % y of point 1
    
    p2x = tmpx(cell_inds == ic & drug_inds == 3); % x of point 2
    p2y = alliqrs(cell_inds == ic & drug_inds == 3)'; % y of point 2
    
    optdir = 1./sqrt(2).*[-1; -1];
    ligdir = [p2x-p1x p2y-p1y];
    score_lengths = diag(ligdir*repmat(optdir,1,size(ligdir,1))); % positive for angles < 90; negative otherwise
    score_angles = pi-acos(score_lengths./(sqrt(ligdir(:,1).^2+ligdir(:,2).^2)));
    
    ctscore = [ctscore sum(score_lengths.*score_angles)];
    
end

figure
bar(ctscore)
set(gca,'XTick',1:length(cell_name),'XTickLabel',cell_name)
ylabel('Diagonal score')


% -------------------------------------------------------------------------
% R-square is no good metric in this case
% close all
% figure
% rsq = [];
% for ic = 1:length(cell_name)
%     rsq = [rsq corr([tmpx(markcond(markcond(:,2) == ic,1)); tmpx(markcond(markcond(:,2) == ic,1)+2)],[alliqrs(markcond(markcond(:,2) == ic,1))'; alliqrs(markcond(markcond(:,2) == ic,1)+2)'])^2];
% end
% bar(rsq)
% set(gca,'XTick',1:length(cell_name),'XTickLabel',cell_name)
% ylabel('R-square')