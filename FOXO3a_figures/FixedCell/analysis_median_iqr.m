%% Immunostaining data of fixed parental cells
cd ..
% cell_name = 'BT20';
% cell_name = 'HS578T';
% cell_name = 'MCF7';
%cell_name = 'MCF10A10142013';
 cell_name = 'MDA231';
% cell_name = 'SKBR3';
% cell_name = 'T47D';
% cell_name = 'HCC1806';

Ligand = {'EGF';'IGF1';'FGF1';'HRG';'HGF';'EPR';'BTC';'NS'}; %row 1-8
Tx = {'No Drug';'AKTi';'MEKi';'Both'}; % col 1-4
obs = {'FOXO3a', 'pERK', 'pAKT'};
timepoints = [0, 15, 30, 60, 90, 120,180,240];

close all

nrows = length(Ligand);
ncols = length(Tx);

% for iobs = 1:length(obs)
for iobs = 1 % Only FOXO3a
    close all
    
    medians = [];
    iqrs = [];
    
    data = load(['./Workspaces/' cell_name '_' obs{iobs}]);
    mydata = getfield(data,['single_' obs{iobs}]);
    
    subplot(nrows,ncols,1)
    text(0,25,0,[cell_name ' - ' obs{iobs}])
    
    posFig = get(gcf,'Position');
    posFig(3) = posFig(3)/2;
    % posFig(4) = posFig(4)*2;
    set(gcf,'Position',posFig)
    set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./10);
    
    for ilig = 1:nrows
        for idrug = 1:ncols
            subplot(nrows,ncols,(ilig-1)*ncols+idrug)
            hold on
            
            for it = 1:length(timepoints)
                data_kern = mydata{it,ilig,idrug}(:,1);
                edges = linspace(min(data_kern),max(data_kern),31);
                plot3(edges,ones(1,length(edges))*it,histc(data_kern,edges))
                view([-5 60])
                
                medians = [medians median(data_kern)];
                iqrs = [iqrs iqr(data_kern)];
            end
            title([Ligand{ilig} ' - ' Tx{idrug}])
            
        end
    end
%     saveFigure(gcf,[cell_name '_' obs{iobs} '_trajectories'],0)
    
    % Analysis of mean / variance
%     close all
    figure

    nrows = length(Ligand);
    ncols = length(Tx);

    posFig = get(gcf,'Position');
    posFig(3) = posFig(3)/2;
    % posFig(4) = posFig(4)*2;
    set(gcf,'Position',posFig)
    set(gcf,'PaperPosition', [0 0 posFig(3) posFig(4)]./10);
    
    subplot(nrows,ncols,1)
    text(0,.3,0,[cell_name ' - ' obs{iobs}])

    for ilig = 1:nrows
        for idrug = 1:ncols
            subplot(nrows,ncols,(ilig-1)*ncols+idrug)

            if ilig == 1 && idrug == 1
                ylabel('iqr')
            end
            if ilig == length(Ligand) && idrug == 1
                xlabel('median')
            end
            hold on

            highlight = [ilig idrug];
            plot(medians,iqrs,'.','Color',[.7 .7 .7])

            if iobs == 1
                % FOXO3a --> quadratic
            else
                % pAKT or pERK --> linear
                axb = polyfit(medians,iqrs,1);
                plot([min(medians) max(medians)],[min(medians) max(medians)]*axb(1) + axb(2),'k--')
            end

            lh = [];
            legstr = {};


            colors = jet(length(timepoints));
            ind = length(timepoints)*((length(Tx)*(highlight(1)-1) + highlight(2)-1))+1;
            ind = ind:(ind+length(timepoints)-1);
            for ih = 1:length(timepoints)
                lh = [lh plot(medians(ind(ih)),iqrs(ind(ih)),'x','Color',colors(ih,:),'LineWidth',3)];
    %             legstr{end+1} = [Ligand{highlight(1)} ' - ' Tx{highlight(2)} '; t = ' num2str(timepoints(ih)) ' min'];
            end
            title([Ligand{highlight(1)} ' - ' Tx{highlight(2)}])

            set(gca,'XLim',[min(medians)-.1 max(medians)+.1],'YLim',[min(iqrs)-.02 max(iqrs)+.02])

        end
    end

    set(gca,'CLim',[0 1])
    colorbar('YTick',linspace(0,1,length(timepoints)),'YTickLabel',timepoints)
    
%     saveFigure(gcf,[cell_name '_' obs{iobs} '_med_iqr'],0)

end


