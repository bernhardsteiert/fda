%% Immunostaining data of fixed parental cells
% cell_name = {'184A1','BT20','HCC1806','HS578T','MCF7','MCF10A','MDA231','SKBR3','T47D'}; % FOXO3a of MDA231 seems not to be normalized
cell_name = {'184A1','BT20','HCC1806','HS578T','MCF7','MCF10A','SKBR3','T47D'};
% cell_name = {'184A1','BT20','HCC1806','HS578T','MCF7','MCF10A','MCF10A12202013','MDA231','SKBR3','T47D'};
% 'MCF10A11292013' == 'MCF10A'

Ligand = {'EGF';'IGF1';'FGF1';'HRG';'HGF';'EPR';'BTC';'NS'}; %row 1-8
Tx = {'No Drug';'AKTi';'MEKi';'Both'}; % col 1-4
obs = {'FOXO3a', 'pERK', 'pAKT'};
timepoints_184a1 = [0, 5, 10, 15, 20, 30, 45, 60, 90, 120, 180, 300, 480];
timepoints_all   = [0,        15,     30,     60, 90, 120, 180, 240];

close all

nColor = length(cell_name);

figure

hold on

colmap = lines(2*nColor);
legh = [];

nrows = length(Ligand);
ncols = length(Tx);

nrowssp = 2;
ncolssp = 4;

% highlight_ligs = [1 5 7];
highlight_ligs = 1:length(Ligand); % this condition is needed for counting to work; others just for plotting

% puls_ligands = nan(length(cell_name),nrows); % Decide for each celltype & ligand combination whether pulsing is possible
puls_ligands = [];
% counter = 0;
for ic = 1:length(cell_name)
% for ic = 1
    subplot(nrowssp,ncolssp,ic)

    if ic == 1
        timepoints = timepoints_184a1;
        tlt80 = 8;
    else
        timepoints = timepoints_all;
        tlt80 = 4;
    end
    
%     for iobs = 1:length(obs)
    for iobs = 1 % Only FOXO3a

        medians = [];
        iqrs = [];

        data = load(['./Workspaces/' cell_name{ic} '_' obs{iobs}]);
        mydata = getfield(data,['single_' obs{iobs}]);

        shiftmedians = [];
        shiftiqrs = [];
        igfhighmedians = [];
        igfhighiqrs = [];

        highlightinds = [];
        highlighttype = [];
        highlightdrug = [];
        for ilig = 1:nrows
            
            for idrug = 1:ncols
                for it = 1:length(timepoints)
                    data = mydata{it,ilig,idrug}(:,1);
                    
                    medians = [medians median(data)];
                    iqrs = [iqrs iqr(data)];
                    
%                     counter = counter +1;
                end
                
%                 if idrug == 1 && ismember(ilig,highlight_ligs)
                if ismember(ilig,highlight_ligs)
                    highlightinds = [highlightinds length(medians)-length(timepoints)+1+tlt80:length(medians)];
                    highlighttype = [highlighttype ones(1,length(timepoints)-tlt80)*ilig];
                    highlightdrug = [highlightdrug ones(1,length(timepoints)-tlt80)*idrug];
                end
                
%                 if idrug == 4 % Both
%                     shiftmedians = [shiftmedians medians(end-length(timepoints)+1:end)]; % all timepoints
%                     shiftiqrs = [shiftiqrs iqrs(end-length(timepoints)+1:end)]; % all timepoints
%                 end
%                 
%                 if ilig == 2 && idrug == 1
%                     igfhighmedians = [igfhighmedians medians(end-length(timepoints)+2:end)]; % neglect first timepoint
%                     igfhighiqrs = [igfhighiqrs iqrs(end-length(timepoints)+2:end)];
%                 end
            end
        end

        % Shift / Rotation based on median (more unbiased to AKTi / IGF)
        [shiftmedians indsorted] = sort(medians);
        shiftiqrs = iqrs(indsorted);
        igfhighmedians = shiftmedians(end-7:end);
        igfhighiqrs = shiftiqrs(end-7:end);
        shiftiqrs = shiftiqrs(1:8);
        shiftmedians = shiftmedians(1:8);
        
        shiftpar = median(shiftmedians);
        shiftiqr = median(shiftiqrs);
        
        medians = medians - shiftpar;
        iqrs = iqrs - shiftiqr;
        
        igfhighmedians = igfhighmedians - shiftpar;
        igfhighiqrs = igfhighiqrs - shiftiqr;
        
        alpha = atan(median(igfhighiqrs)/median(igfhighmedians));
        Rmat = [cos(alpha)  sin(alpha); ...
                -sin(alpha) cos(alpha)];
        tmp = Rmat*[medians; iqrs];
        medians = tmp(1,:);
        iqrs = tmp(2,:);
        
        igfhighmedians = sort(medians);
        igfhighmedians = igfhighmedians(end-7:end);
        
        medians = medians./median(igfhighmedians);
        
        
%         legh = [legh plot(medians,iqrs,'o','Color',colmap(ic,:))];
        plot(medians,iqrs,'ko')
        title(cell_name{ic})
        hold on
        x = linspace(0,1,201);
        iqrcut = iqrs > max(iqrs)*.2;
        b = lsqnonlin(@(b) b*medians(iqrcut).^2-b*medians(iqrcut) - iqrs(iqrcut),-.5,-Inf,0,optimset('Display','off'));
        y = b*x.^2-b*x;
        plot(x,y,'Color',colmap(ic,:))
        xgate = .2; % must be between 0 and .5
        [tmp ind] = min(abs(x-xgate));
        center = .5;
        plot([center x(ind)],[0 y(ind)],'k-')
        plot([center 2*center-x(ind)],[0 y(ind)],'k-')
        
        % Detecting conditions that lie within 'cone'
        slope = abs(y(ind)./(center-x(ind)));
        indpuls = abs(iqrs./(center-medians)) > slope & iqrs > 0;
%         legh = [legh plot(medians(indpuls),iqrs(indpuls),'o','MarkerEdgeColor','k','MarkerFaceColor',colmap(ic,:))];
        indpuls = find(indpuls);
        
        % Highlighting late behavior (t > 80min) for defined ligands
        for itype = 1:length(highlight_ligs)
            for idrug = 1:ncols
                myinds = highlighttype == highlight_ligs(itype) & highlightdrug == idrug;
                if idrug == 1
                    legh = [legh plot(medians(highlightinds(myinds)),iqrs(highlightinds(myinds)),'o','MarkerFaceColor',colmap(itype,:))];
                end

                members = ismember(indpuls,highlightinds(myinds));
    %             puls_ligands(ic,itype) = sum(members) > (sum(highlighttype == highlight_ligs(itype))*.5);
                puls_ligands = [puls_ligands; sum(members) > (sum(myinds)*.8)];
            end
        end
        
        legend(legh,Ligand{highlight_ligs})

    end
    
    set(gca,'XLim',[-.2 1.2],'YLim',[-.05 .35])

end
% legend(legh,cell_name)
