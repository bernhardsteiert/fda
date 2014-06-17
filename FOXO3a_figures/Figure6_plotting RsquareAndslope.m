close all;clear all;
load rawdata
load('mycolormap','mycmap')
testParam = 'MEKidrivenFOXO3aiqrchange';%'AKTidrivenFOXO3amediancha';%
AllCellLines = unique(cellstr(rawdata(:,'Cell')));
plotCells = {'HCC1806';'MCF10A'};
mycolor = hsv(length(plotCells));
lind = 1;
plotInd = 1;
myLegend = [];
for i=1:length(AllCellLines)

    cellline = AllCellLines{i};
    plotCriteria = rawdata.Cell==cellline ;
    myX = rawdata(plotCriteria,'medianpAKT');
    myY = rawdata(plotCriteria,'medianFOXO3a');
    
    myColor = rawdata(plotCriteria,testParam);
    figure(i),h = scatter(myX,myY,60,myColor,'fill','o');
    xlabel('Median pAKT');
    ylabel('Median FoxO3a');
    title(cellline);
    set(h,'MarkerEdgeColor','none');
    set(gca,'CLim',[-2 2],'XTick',0:0.25:1,'YTick',0:0.25:1);
    set(i,'Colormap',mycmap);
    xlim([-0.1 1.1]);
    ylim([-0.1 1.1]);
    clear Xn Yn p f x1 x2 X z IX b x1fit x2fit YFIT X1FIT X2FIT;
    
    [Xn,IX] = sort(double(myX));
    Yn = double(myY);
    Yn = Yn(IX);
    p = polyfit(Xn,Yn,2);
    f = polyval(p,Xn);
    hold on;plot(Xn,f,'-k'); hold off;
    
    x1 = Xn;
    x2 = Yn;
    X = [ones(size(x1)) x1 x2 x1.*x2];
    z = double(myColor);
    z = z(IX);
    [b,bint,r,rint,stats] = regress(z,X);
    Rsquare1(i) = stats(1);
    Fstat1(i) = stats(2);
    PVal1(i) = stats(3);
    B2(i) = b(2);
    B3(i) = b(3);
    B4(i) = b(4);
    
    figure(i+10);
    scatter3(x1,x2,z,'filled'); hold on;

    x1fit = min(x1):(max(x1)-min(x1))/10:max(x1);
    x2fit = min(x2):(max(x2)-min(x2))/10:max(x2);
    [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
    YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
    mesh(X1FIT,X2FIT,YFIT)
    xlabel('Median pAKT');
    ylabel('Median FoxO3a');
    zlabel(' log_2(FoxO3a IQR change)');
    title([cellline ', R^2=' num2str(stats(1))]);
    view(39,12)
    
    
    TPs = unique(rawdata.Time);
    storeX = [];
    storeY = [];
    Rsquare = [];
    Fstat = [];
    PVal = [];
    Slope = [];
    
    for t = 1:length(TPs)
        plotCriteria = rawdata.Cell==cellline & rawdata.Time == TPs(t) %& rawdata.MEKidrivenFOXO3aiqrchange < 0 & rawdata.MEKidrivenpAKTmedianchang > 0;
        myX = double(rawdata(plotCriteria,'MEKidrivenpAKTmedianchang'));
        myY = double(rawdata(plotCriteria,testParam));
        [b,bint,r,rint,stats] = regress(myY,[ones(size(myX)),myX]);
        Rsquare(t) = stats(1);
        Fstat(t) = stats(2);
        PVal(t) = stats(3);
        Slope(t) = b(2);
        %if   Rsquare(t) > 0.3 
            storeX = [storeX;myX];
            storeY = [storeY;myY];
        %end
    end
    
    if ~isempty(storeX)
        [b,bint,r,rint,stats] = regress(storeY,[ones(size(storeX)),storeX]);
        
        for cInd = 1:length(plotCells)
            if strcmp(cellline,plotCells{cInd})
                figure(22);hold on;plot(storeX,storeY,'.','color',mycolor(plotInd,:),'MarkerSize',20);hold off;
                myLegend{plotInd} = cellline;plotInd=plotInd+1;
            end
        end
        
        Rsquare2(i) = stats(1);
        Fstat2(i) = stats(2);
        PVal2(i) = stats(3);
        Slope2(i) = b(2);
        %YFIT = b(1) + b(2)*storeX ;
        %hold on; plot(storeX,YFIT,'-','color',mycolor(lind,:));
        
        lind = lind+1;
    else
        Rsquare2(i) = 0;
        Fstat2(i) = 0;
        PVal2(i) = 0;
        Slope2(i) = 0;
    end
    clear storeX storeX stats
end

figure(22);legend(myLegend);
xlabel('pAKT Median change');
ylabel('log_2(FoxO3a IQR change)');
Fstat1(isnan(Fstat1)) = 0;
[sortedFstat1 IDX] = sort(Fstat1);

figure(),
subplot(6,1,1);bar(Rsquare1(IDX)); set(gca,'XTickLabel',AllCellLines(IDX));title('R^2');
subplot(6,1,2);bar(sortedFstat1); set(gca,'XTickLabel',AllCellLines(IDX));title('F-stat');
subplot(6,1,3);bar(PVal1(IDX)); set(gca,'XTickLabel',AllCellLines(IDX));title('P-value');
subplot(6,1,4);bar(B2(IDX)); set(gca,'XTickLabel',AllCellLines(IDX));title('b*X');
subplot(6,1,5);bar(B3(IDX)); set(gca,'XTickLabel',AllCellLines(IDX));title('b*Y');
subplot(6,1,6);bar(B4(IDX)); set(gca,'XTickLabel',AllCellLines(IDX));title('b*XY');
clear IDX;

%Slope2(Rsquare2 < 0.3) = 0;
%Fstat2(Rsquare2 < 0.3) = 0;
%PVal2(Rsquare2 < 0.3) = 0;
%Rsquare2(Rsquare2 < 0.3) = 0;

[sortedSlope2 IDX] = sort(-Slope2);
figure(),
subplot(4,1,1);bar(Rsquare2(IDX)); set(gca,'XTickLabel',AllCellLines(IDX));title('R^2');
subplot(4,1,2);bar(Fstat2(IDX)); set(gca,'XTickLabel',AllCellLines(IDX));title('F-stat');
subplot(4,1,3);bar(PVal2(IDX)); set(gca,'XTickLabel',AllCellLines(IDX));title('P-value');
subplot(4,1,4);bar(Slope2(IDX)); set(gca,'XTickLabel',AllCellLines(IDX));title('Slope');

figure(),
bar(Rsquare2(IDX)); set(gca,'XTickLabel',AllCellLines(IDX));title('R^2');
figure(),
bar(-Slope2(IDX)); set(gca,'XTickLabel',AllCellLines(IDX));title('Slope');