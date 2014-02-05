% Plot EGF early vs. late time-points
close all
clear all

log_trafo = 1;

load('./Workspaces/site_64')
c_signal = log10(intensity);

time_range = [50 510];

xpanpos = -.15;
ypanpos = 1.05;
fontsize = 16;

[tmp range_ind_min] = min(abs(timestamp - time_range(1)));
[tmp range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

f1 = figure;

xfac = 1.5;
yfac = 1;
fontsizef = 10;

setFigure(f1,xfac,yfac,fontsizef)

nrows = 2;
ncols = 3;
% subplot(nrows,ncols,[1 2 3 6 7 8 11 12 13])
s1 = subplot(nrows,ncols,1);
s1Pos = get(s1,'Position');

% Exclude Inhibitor data
isite = 64;
s = siteprop(isite);
site_lig_ind = s.lig_index;
site_lig_name = s.lig_name;
site_lig_dose = s.lig_dose;
site_inh_name = s.inh_name;
site_inh_dose = s.inh_dose;

first_n = 7; % Plot first_n traces colored (2 times)
set(f1,'DefaultAxesColorOrder',lines(first_n))

tmp_puls_strengths = edge_snr_score_pw_distdur(isite);
[tmp ind_tmp_pul_str] = sort(tmp_puls_strengths);
ind_isite = [ind_tmp_pul_str(end:-1:end-first_n+1) ind_tmp_pul_str(round(linspace(1,(length(ind_tmp_pul_str)-first_n),first_n)))];

plot(repmat(timestamp(range_ind),1,size(c_signal,2)),c_signal(range_ind,:),'g','color',[0.7 0.7 0.7])
hold on
plot(repmat(timestamp(range_ind),1,2*first_n),c_signal(range_ind,ind_isite))
plot(timestamp(range_ind),nanmean(c_signal(range_ind,:),2),'color','k','LineWidth',2)
title([site_lig_name num2str(site_lig_dose) ' ng/ml'])

xlabel('time [min]')
ylabel('log_{10} FOXO3a [Cyt/Nuc]');

ylim = [-1 1]*.04;
if ~log_trafo
    ylim = 10.^ylim;
end
plot([200 200],ylim,'k--')
text(125,.033,{'non-stationary','(deterministic)'},'HorizontalAlignment','center')
text(355,.033,{'stationary','(stochastic)'},'HorizontalAlignment','center')

set(gca,'XLim',time_range,'YLim',ylim)
set(gca,'XTick',50:50:500)

text(xpanpos, ypanpos,'A','units', 'normalized','FontSize',fontsize)

% Figure 2B: Eigenfunctions (new - rotated)
subplot(nrows,ncols,2)
hold on

load('harm_basis_fPCA.mat') % Contains only harm_basis_fPCA from all data-sets

% time_range = getbasisrange(harm_basis_fPCA);
time_range = [50.75 197.91];
times_fine = linspace(time_range(1),time_range(2),501);

basis_eval = eval_basis(harm_basis_fPCA,times_fine);

ylim = [0 1];
xlim = [50 510];
set(gca,'XLim',[50 510],'YLim',ylim)
plot([198 198],ylim,'k')
plot([202 202],ylim,'k')
plot(xlim,[1/3 1/3],'k')
plot(xlim,[2/3 2/3],'k')
set(gca,'XTick',50:50:500)
set(gca,'YTick',[1/6 3/6 5/6],'YTickLabel',[3 2 1],'TickLength',[0 0])
box on
xlabel('time [min]')
ylabel('Harmonic')

range_eval = 2.2*max(max(abs(basis_eval)));

for iplot = 1:size(basis_eval,2)
%     subplot(nrows,ncols,ncols*(iplot-1)+ncols-1)
    
    tmpplot = basis_eval(:,iplot)./(range_eval*3)+((3-iplot)*2+1)/6;
    
    plot(times_fine,tmpplot,'k')   
    plot(time_range,[((iplot-1)*2+1)/6 ((iplot-1)*2+1)/6],'k--')
    
end

load('harm_basis.mat') % Contains only harm_basis from all data-sets

% time_range = getbasisrange(harm_basis);
time_range = [202 510];
times_fine = linspace(time_range(1),time_range(2),501);
basis_eval = eval_basis(harm_basis,times_fine);

% range_eval = 2.8*max(max(abs(basis_eval)));

for iplot = 1:size(basis_eval,2)
%     subplot(nrows,ncols,ncols*(iplot-1)+ncols-1)
    
    tmpplot = basis_eval(:,iplot)./(range_eval*3)+((3-iplot)*2+1)/6;
    
    plot(times_fine,tmpplot,'k')   
    plot(time_range,[((iplot-1)*2+1)/6 ((iplot-1)*2+1)/6],'k--')
    
end

text(xpanpos, ypanpos,'B','units', 'normalized','FontSize',fontsize)


s3 = subplot(nrows,ncols,3);
hold on
sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:70]; % Without FGF
pcs = [2 3];
possible_doses = [0 2.5 5 10 20 50 100];

el_area = nan(6,7); % 6 ligands; 7 doses
var_pul = nan(6,7);

highdoses = [];
resort = [4 1 nan 2 3 6 5]; % Relative to platemap
for isite = sites_all
    scores = fPCA(isite);
    plot(scores(2,:),scores(3,:),'.','Color',[.7 .7 .7])
    
    sprop = siteprop(isite);
    if sprop.lig_dose == 100
        highdoses = [highdoses isite];
    end
    
    el_area(resort(sprop.lig_index),sprop.lig_dose == possible_doses) = ellipsis_area(scores(2,:),scores(3,:));
    var_pul(resort(sprop.lig_index),sprop.lig_dose == possible_doses) = iqr(edge_snr_score_pw_distdur(isite));
end

resort = [2 3 4 1 6 5];
% resort = [6:-1:1];
highdoses = highdoses(resort);

color_ind = 1;
colmap = [linspace(0,1,length(highdoses)+1)' ones(length(highdoses)+1,1) ones(length(highdoses)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
legstr = cell(length(highdoses),1);
for isite = highdoses([6 2 3 4 5 1])
    sprop = siteprop(isite);
    legstr{isite == highdoses} = sprop.lig_name(1:3);
    
    scores = fPCA(isite);
    plot(scores(2,:),scores(3,:),'o','MarkerSize',3,'MarkerFaceColor',colmap(isite == highdoses,:),'MarkerEdgeColor','none')
    plotEllipsis(scores(2,:),scores(3,:),colmap(isite == highdoses,:),.5);
end

xlim = [-.16 .24];
ylim = [-.1 .14];
set(gca,'XLim',xlim,'YLim',ylim)

% axisEqual(get(gcf,'Position'))

ylabel(['Score PC' num2str(pcs(2))])
% arrow([-.13 -.1],[.2 -.1],'Width',.5,'Length',7)
% set(gca,'YTick',-.1:.1:.2)
xlabel(['Score PC' num2str(pcs(1))])
% arrow([-.13 .06],[.02 .15],'Width',.5,'Length',7)

set(gca,'CLim',[0 1])
subplotpos = get(s3,'Position');
colormap(colmap)
colorbar('YTick',linspace(1./(2*length(highdoses)),1-1./(2*length(highdoses)),length(highdoses)),'YTickLabel',legstr,'TickLength', [0 0],'Position',[subplotpos(1)+subplotpos(3) subplotpos(2) .01 subplotpos(4)],'units','normalized') % Vertical colorbar
box on

text(xpanpos, ypanpos,'C','units', 'normalized','FontSize',fontsize)

s4 = subplot(nrows,ncols,4);
sites = [4 10 17 37 44 57 64];
input_names = {'EGF','No Lig','IGF','HRG','HGF','EPR','BTC'};

dists = [];
pca_scores = [];
celltypeharm = [];

dists_mat = [];

for isite = sites
    radial_dists = edge_snr_score_pw_distdur(isite);
    
    dists = [dists radial_dists];
    pca_scores = [pca_scores fPCA(isite)];
    
    celltypeharm = [celltypeharm ones(size(radial_dists))*isite];
    
    dists_mat = padconcatenation(dists_mat,radial_dists,1);

end

resort = [3 2 6 4 1 5 7];
boxplot(dists_mat(resort,:)')

text(xpanpos, ypanpos,'D','units', 'normalized','FontSize',fontsize)

s4Pos = get(s4,'Position');
set(s4,'Position',[s1Pos(1) s4Pos(2) s1Pos(3:4)])

ylabel('Pulsatory strength')
set(gca,'XTick',1:7,'XTickLabel',input_names(resort))

s5 = subplot(nrows,ncols,5);
colmap2 = lines(3);
colmap2(1,:) = [0 0 0];
color_ind = 1;
leghand = [];
hold on
% for ipc = 1:3
%     mycolor = colmap2(color_ind,:);
%     color_ind = color_ind + 1;
%     for isite = sites
%         % Colored
%         scores = fPCA(isite);
%         if find(isite==sites) == 1
%             leghand = [leghand plot(dists(celltypeharm == isite),scores(ipc,:),'o','MarkerFaceColor',mycolor,'MarkerEdgeColor','none','MarkerSize',3)];
%         else
%             plot(dists(celltypeharm == isite),scores(ipc,:),'o','MarkerFaceColor',mycolor,'MarkerEdgeColor','none','MarkerSize',3)
%         end
%     end
% end
% legend(leghand,{'PC1','PC2','PC3'})
% set(gca,'XLim',[0 1.6],'XTick',0:.2:1.6)
% set(gca,'YLim',[-.25 .25])

[bcoef,bint,r,rint,stats] = regress(dists',[pca_scores' ones(size(pca_scores,2),1)]);
plot(dists,[pca_scores' ones(size(pca_scores,2),1)]*bcoef,'o','MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',3)
set(gca,'XLim',[0 1.6],'XTick',0:.2:1.6)
set(gca,'YLim',[-.1 1.1])
plot([0 1],[0 1],'k','LineWidth',2)

text(1.1,.9,sprintf('$R^2=%1.3f$',stats(1)),'Interpreter','latex')
xlabel('Pulsatory strength')
ylabel('$\sum_i b_i s_i + c$','Interpreter','latex')
box on

s4Pos = get(s4,'Position');
s5Pos = get(s5,'Position');
set(s5,'Position',[s5Pos(1) s4Pos(2) s5Pos(3:4)])

text(xpanpos, ypanpos,'E','units', 'normalized','FontSize',fontsize)

s6 = subplot(nrows,ncols,6);
hold on
% legh = [];
xMinMax = [min(min(var_pul))-.1 max(max(var_pul))+.1];
yMinMax = [min(min(el_area))-.002 max(max(el_area))+.002];
set(gca,'XLim',xMinMax)
set(gca,'YLim',yMinMax)

s6Pos = get(s6,'Position');
set(s6,'Position',[s6Pos(1) s4Pos(2) s6Pos(3:4)])

axPos = get(gca,'Position');
% for ilig = 1:size(el_area,1)
%     for idose = 1:size(el_area,2)-1
%         xAnnotation = axPos(1) + ((var_pul(ilig,idose:idose+1) - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
%         yAnnotation = axPos(2) + ((el_area(ilig,idose:idose+1) - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);
% 
%         annotation('arrow',xAnnotation,yAnnotation,'Color',colmap(ilig,:))
%     end
% end
% % legend(legh,legstr)

for ilig = 1:size(el_area,1)
    xAnnotation = axPos(1) + ((var_pul(ilig,end) - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
    yAnnotation = axPos(2) + ((el_area(ilig,end) - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);
    
    plot(var_pul(ilig,end),el_area(ilig,end),'o','MarkerFaceColor',colmap(ilig,:))
    annotation('textbox',[xAnnotation yAnnotation 0 0],'String',legstr{ilig},'LineStyle','none','HorizontalAlignment','center')
end
xAnnotation = axPos(1) + ((mean(var_pul(:,1)) - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
yAnnotation = axPos(2) + ((mean(el_area(:,1)) - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);

plot(mean(var_pul(:,1)),mean(el_area(:,1)),'o','MarkerFaceColor','k')
offset_an = .06;
annotation('textbox',[xAnnotation-offset_an/2 yAnnotation offset_an 0],'String','No Lig','LineStyle','none','HorizontalAlignment','center')

xlabel('IQR Pulsatory strength')
ylabel('Area ellipse')
box on

text(xpanpos, ypanpos,'F','units', 'normalized','FontSize',fontsize)