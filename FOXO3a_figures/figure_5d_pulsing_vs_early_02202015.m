% Figure 5d: Pulsing vs. early PC2
% New dataset: 02202015.

load('Workspaces/scores_02202015')

pc = 2;

addpath('./Functions/')

extension = '02-20-2015_cleaned';
sites = 1:60;

puls_thres = .55;

meki_doses = [.125/3^0 .125/3^1 .125/3^2 .125/3^3 .125/3^4 0];
akti_doses = [1/1.5^0 1/1.5^1 1/1.5^2 1/1.5^3 1/1.5^4 1/1.5^5 1/1.5^6 1/1.5^7 1/1.5^8 1/1.5^9 1/1.5^10 0];

egfmeki_pulsing = [];
egfmeki_early = [];
egfmeki_akti_dose = [];
egfmeki_meki_dose = [];
reorder = nan(1,length(sites));

for isite = sites

    s = siteprop(isite,extension);

    idrug1 = s.drug1_dose;
    idrug2 = s.drug2_dose;

    egfmeki_pulsing = [egfmeki_pulsing sum(dists(celltype == isite) > puls_thres) / sum(celltype == isite)];
    egfmeki_early = [egfmeki_early nanmean(scores_all(pc,celltype == isite)-scores_all(1,celltype == isite))];
    egfmeki_akti_dose = [egfmeki_akti_dose idrug1];
    egfmeki_meki_dose = [egfmeki_meki_dose idrug2];
    reorder(isite) = subplotpos(isite,12);

end

egfmeki_early_mat = nan(12,5);
egfmeki_early_mat(:) = egfmeki_early(reorder);
egfmeki_early_mat = egfmeki_early_mat';

egfmeki_pulsing_mat = nan(12,5);
egfmeki_pulsing_mat(:) = egfmeki_pulsing(reorder);
egfmeki_pulsing_mat = egfmeki_pulsing_mat';

figure
hold on
title('Fixed MEKi for each color, larger marker = less AKTi')

for i = sites
    markersize = 2 + find(egfmeki_akti_dose(i) == akti_doses) * 2.5;
    col2 = .7;
    col1 = (length(meki_doses)-find(egfmeki_meki_dose(i) == meki_doses)+1) / length(meki_doses);
    col3 = 1;
    doplot = 1;
    if doplot
        col = hsv2rgb([col1 col2 col3]);
        plot(egfmeki_early(i),egfmeki_pulsing(i),'o','MarkerSize',markersize,'MarkerFaceColor',col,'MarkerEdgeColor','none')
    end
end

legh = [];
legstr = {};
for i = 1:size(egfmeki_early_mat,1)
    col2 = .7;
    col1 = (length(meki_doses)-i+1) / length(meki_doses);
    col3 = 1;
    col = hsv2rgb([col1 col2 col3]);
    legh = [legh plot(egfmeki_early_mat(i,:),egfmeki_pulsing_mat(i,:),'-.','Color',col)];
    legstr{end+1} = sprintf('MEKi = %g',1000*meki_doses(i));
end
legend(legh,legstr)

xlabel(sprintf('Early PC%i',pc))
ylabel('Fraction of pulsing cells')

figure
hold on
title('Fixed AKTi for each color; larger marker = less MEKi')

for i = sites
    markersize = 4 + find(egfmeki_meki_dose(i) == meki_doses) * 4;
    col2 = .7;
    col1 = (length(akti_doses)-find(egfmeki_akti_dose(i) == akti_doses)+1) / length(akti_doses);
    col3 = 1;
    doplot = 1;
    if doplot
        col = hsv2rgb([col1 col2 col3]);
        plot(egfmeki_early(i),egfmeki_pulsing(i),'o','MarkerSize',markersize,'MarkerFaceColor',col,'MarkerEdgeColor','none')
    end
end

legh = [];
legstr = {};
for i = 1:size(egfmeki_early_mat,2)
    col2 = .7;
    col1 = (length(akti_doses)-i+1) / length(akti_doses);
    col3 = 1;
    col = hsv2rgb([col1 col2 col3]);
    legh = [legh plot(egfmeki_early_mat(:,i),egfmeki_pulsing_mat(:,i),'-.','Color',col)];
    legstr{end+1} = sprintf('AKTi = %g',1000*akti_doses(i));
end
legend(legh,legstr)

xlabel(sprintf('Early PC%i',pc))
ylabel('Fraction of pulsing cells')

