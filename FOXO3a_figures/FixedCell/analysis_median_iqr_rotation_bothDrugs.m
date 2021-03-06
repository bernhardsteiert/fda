%% Immunostaining data of fixed parental cells
cell_name = {'184A1','MCF10A11292013','SKBR3','HCC1806','HS578T','MDA231','BT20','MCF7','T47D'};
%cell_name = {'184A1','MCF10A11292013','SKBR3','HCC1806','HS578T','MDA231','BT20','MCF7','T47D'};

% 'MCF10A11292013' == 'MCF10A'

Ligand = {'EGF';'IGF1';'FGF1';'HRG';'HGF';'EPR';'BTC';'NS'}; %row 1-8
Tx = {'No Drug';'AKTi';'MEKi';'Both'}; % col 1-4
obs = {'FOXO3a', 'pERK', 'pAKT'};
timepoints_184a1 = [0, 5, 10, 15, 20, 30, 45, 60, 90, 120, 180, 300, 480];
timepoints_all   = [0,        15,     30,     60, 90, 120, 180, 240];

highlight_ligs = [1]; % EGF; FGF1; BTC; NS
%highlight_ligs = [2 4 5]; % IGF1 HRG NS
testDrug = 3;
close all

figure
set(gcf,'Position',[100 100 984 638]);
hold on

colmap = [hsv(length(highlight_ligs))];

nrows = length(Ligand);
ncols = length(Tx);

nrowssp = 3;
ncolssp = 3;

markers = {'o','x','^'};

for ic = 1:length(cell_name)

    subplot(nrowssp,ncolssp,ic)

    if ic == 1
        timepoints = timepoints_184a1;
    else
        timepoints = timepoints_all;
    end
    
    for iobs = 1 % Only FOXO3a

        medians = [];
        iqrs = [];
        currentLoc = pwd;
        cd ..;
        data = load(['./Workspaces/' cell_name{ic} '_' obs{iobs}]);
        cd(currentLoc);
        mydata = getfield(data,['single_' obs{iobs}]);

        shiftmedians = [];
        shiftiqrs = [];
        igfhighmedians = [];
        igfhighiqrs = [];

        allind = [];
        alltype = [];
        alldrug = [];
        alltime = [];
        for ilig = 1:nrows
            
            for idrug = 1:ncols
                for it = 1:length(timepoints)
                    data = mydata{it,ilig,idrug}(:,1);
                    
                    medians = [medians median(data)];
                    iqrs = [iqrs iqr(data)];
                    
                    allind = [allind length(allind)+1];
                    alltype = [alltype ilig];
                    alldrug = [alldrug idrug];
                    alltime = [alltime timepoints(it)];
%                     counter = counter +1;
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
        
        angle = atan(median(igfhighiqrs)/median(igfhighmedians));
        Rmat = [cos(angle)  sin(angle); ...
                -sin(angle) cos(angle)];
        tmp = Rmat*[medians; iqrs];
        medians = tmp(1,:);
        iqrs = tmp(2,:);
        
        igfhighmedians = sort(medians);
        igfhighmedians = igfhighmedians(end-7:end);
        
        medians = medians./median(igfhighmedians);
        
        
%         legh = [legh plot(medians,iqrs,'o','Color',colmap(ic,:))];
%         plot(medians,iqrs,'ko')
        x = linspace(0,1,201);
        iqrcut = iqrs > max(iqrs)*.2;
        b = lsqnonlin(@(b) b*medians(iqrcut).^2-b*medians(iqrcut) - iqrs(iqrcut),-.5,-Inf,0,optimset('Display','off'));
        y = b*x.^2-b*x;
%         plot(x,y,'Color',colmap(ic,:))
        plot(x,y,'Color','k')
        title(cell_name{ic})
        hold on
        % Plot cone
        xgate = .2; % must be between 0 and .5
        [tmp ind] = min(abs(x-xgate));
        center = .5;
%         plot([center x(ind)],[0 y(ind)],'k-')
%         plot([center 2*center-x(ind)],[0 y(ind)],'k-')
        
        % Detecting conditions that lie within 'cone'
        slope = abs(y(ind)./(center-x(ind)));
        indpuls = abs(iqrs./(center-medians)) > slope & iqrs > 0;
%         legh = [legh plot(medians(indpuls),iqrs(indpuls),'o','MarkerEdgeColor','k','MarkerFaceColor',colmap(ic,:))];
        indpuls = find(indpuls);
        
        
        % Highlighting late behavior (t > 80min) for defined ligands
        legh = [];
        legstr = {};
        plottedinds = [];
        for itype = 1:length(highlight_ligs)
            idrug = 1;
            myinds = allind(alltype == highlight_ligs(itype) & alldrug == idrug & alltime > 0);
            legh = [legh plot(medians(myinds),iqrs(myinds),markers{1},'MarkerFaceColor','none','MarkerEdgeColor',colmap(itype,:))];
            legstr{end+1} = Ligand{highlight_ligs(itype)};
            plottedinds = [plottedinds myinds];

            idrug = testDrug;
            myinds = allind(alltype == highlight_ligs(itype) & alldrug == idrug & alltime > 0);
            legh = [legh plot(medians(myinds),iqrs(myinds),markers{1},'MarkerFaceColor',colmap(itype,:),'MarkerEdgeColor',colmap(itype,:))];
            legstr{end+1} = [Ligand{highlight_ligs(itype)} '+' Tx{idrug}];
            plottedinds = [plottedinds myinds];
        end
        
        
        
        
%         plottedinds = [plottedinds myinds];
%         myinds = ~ismember(1:length(medians),plottedinds);
%         plot(medians(myinds),iqrs(myinds),'x','MarkerFaceColor','none','MarkerEdgeColor',[.75 .75 .75]);
%         uistack(legh,'top')
% 
        if ic == length(cell_name)
            legend(legh,legstr)
        end

    end
    
    set(gca,'XLim',[-.2 1.2],'YLim',[-.05 .35])
    
end



% for AKT
figure('Name','pAKT');
for ic = 1:length(cell_name)

    subplot(nrowssp,ncolssp,ic)

    if ic == 1
        timepoints = timepoints_184a1;
    else
        timepoints = timepoints_all;
    end
    
    for iobs = 3 % Only pAKT
        title(cell_name{ic});
        medians = [];
        iqrs = [];
        currentLoc = pwd;
        cd ..;
        data = load(['./Workspaces/' cell_name{ic} '_' obs{iobs}]);
        cd(currentLoc);
        mydata = getfield(data,['single_' obs{iobs}]);

        allind = [];
        alltype = [];
        alldrug = [];
        alltime = [];
        for ilig = 1:nrows
            
            for idrug = 1:ncols
                for it = 1:length(timepoints)
                    data = mydata{it,ilig,idrug}(:,1);
                    
                    medians = [medians median(data)];
                    iqrs = [iqrs iqr(data)];
                    
                    allind = [allind length(allind)+1];
                    alltype = [alltype ilig];
                    alldrug = [alldrug idrug];
                    alltime = [alltime timepoints(it)];
                end
            end
        end
        
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
        
        angle = 0;
        Rmat = [cos(angle)  sin(angle); ...
                -sin(angle) cos(angle)];
        tmp = Rmat*[medians; iqrs];
        medians = tmp(1,:);
        iqrs = tmp(2,:);
        
        igfhighmedians = sort(medians);
        igfhighmedians = igfhighmedians(end-7:end);
        
        medians = medians./median(igfhighmedians);
        

        legh = [];
        legstr = {};
        plottedinds = [];
        hold on;
        for itype = 1:length(highlight_ligs)
             
            idrug = testDrug;
            myinds = allind(alltype == highlight_ligs(itype) & alldrug == idrug & alltime > 0);
            legh = [legh plot(medians(myinds),iqrs(myinds),markers{1},'MarkerFaceColor',colmap(itype,:),'MarkerEdgeColor',colmap(itype,:))];
            legstr{end+1} = [Ligand{highlight_ligs(itype)} '+' Tx{idrug}];
            plottedinds = [plottedinds myinds];                  
            
            idrug = 1;
            myinds = allind(alltype == highlight_ligs(itype) & alldrug == idrug & alltime > 0);
            legh = [legh plot(medians(myinds),iqrs(myinds),markers{1},'MarkerFaceColor','none','MarkerEdgeColor',colmap(itype,:))];
            legstr{end+1} = Ligand{highlight_ligs(itype)};
            plottedinds = [plottedinds myinds];
     
        end
        if ic == length(cell_name)
            legend(legh,legstr)
        end
        set(gca,'XLim',[-.1 1.1],'YLim',[-.05 .3])
    end
end


% for ERK
figure('Name','pERK');
for ic = 1:length(cell_name)

    subplot(nrowssp,ncolssp,ic)

    if ic == 1
        timepoints = timepoints_184a1;
    else
        timepoints = timepoints_all;
    end
    
    for iobs = 2 % Only ERK
        title(cell_name{ic});
        medians = [];
        iqrs = [];
        currentLoc = pwd;
        cd ..;
        data = load(['./Workspaces/' cell_name{ic} '_' obs{iobs}]);
        cd(currentLoc);
        mydata = getfield(data,['single_' obs{iobs}]);

        allind = [];
        alltype = [];
        alldrug = [];
        alltime = [];
        for ilig = 1:nrows
            
            for idrug = 1:ncols
                for it = 1:length(timepoints)
                    data = mydata{it,ilig,idrug}(:,1);
                    
                    medians = [medians median(data)];
                    iqrs = [iqrs iqr(data)];
                    
                    allind = [allind length(allind)+1];
                    alltype = [alltype ilig];
                    alldrug = [alldrug idrug];
                    alltime = [alltime timepoints(it)];
                end
            end
        end
        
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
        
        angle = 0;
        Rmat = [cos(angle)  sin(angle); ...
                -sin(angle) cos(angle)];
        tmp = Rmat*[medians; iqrs];
        medians = tmp(1,:);
        iqrs = tmp(2,:);
        
        igfhighmedians = sort(medians);
        igfhighmedians = igfhighmedians(end-7:end);
        
        medians = medians./median(igfhighmedians);
        

        legh = [];
        legstr = {};
        plottedinds = [];
        hold on;
        for itype = 1:length(highlight_ligs)
             
            idrug = testDrug;
            myinds = allind(alltype == highlight_ligs(itype) & alldrug == idrug & alltime > 0);
            legh = [legh plot(medians(myinds),iqrs(myinds),markers{1},'MarkerFaceColor',colmap(itype,:),'MarkerEdgeColor',colmap(itype,:))];
            legstr{end+1} = [Ligand{highlight_ligs(itype)} '+' Tx{idrug}];
            plottedinds = [plottedinds myinds];                  
            
            idrug = 1;
            myinds = allind(alltype == highlight_ligs(itype) & alldrug == idrug & alltime > 0);
            legh = [legh plot(medians(myinds),iqrs(myinds),markers{1},'MarkerFaceColor','none','MarkerEdgeColor',colmap(itype,:))];
            legstr{end+1} = Ligand{highlight_ligs(itype)};
            plottedinds = [plottedinds myinds];
     
        end
        if ic == length(cell_name)
            legend(legh,legstr)
        end
        set(gca,'XLim',[-.1 1.1],'YLim',[-.05 .3])
    end
end

