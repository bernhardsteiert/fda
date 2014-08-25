% Figure supplementary: %variance explained vs. #basis functions
addpath('./Functions/')

load('./Workspaces/harm_basis_130722_corrected_retracked_all_cleaned_late_newBTC')

deltat = 120;

time_range = getbasisrange(harm_basis);
timestamp = linspace(time_range(1),time_range(2),501);
harm_eval = eval_basis(harm_basis,timestamp);
timestamp = timestamp - deltat;

f1 = figure;

xfac = .6;
yfac = 1;
fontsize = 16;

setFigure(f1,xfac,yfac,fontsize)

for i = 1:size(harm_eval,2)
    subplot(3,1,i)

    plot(timestamp,harm_eval(:,i),'k')
    hold on
    plot(timestamp,timestamp*0,'k--')

    xlabel('time [min]')
    ylabel(sprintf('Harmonic %i',i))
    
    set(gca,'XLim',time_range-deltat,'YLim',[-.065 .065])
end