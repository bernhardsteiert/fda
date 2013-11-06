function saveFigure(h, name)

savePath =  ['./Figures/' datestr(now, 30) '_' name];

axis_handle = get(gcf, 'CurrentAxes');
xlabel(axis_handle, 'time [min]')
ylabel(axis_handle, 'FOXO3a Ratio [C/N]')

saveas(h, savePath, 'fig');
print('-depsc', savePath);
% print('-dpng', savePath);
system(['export LD_LIBRARY_PATH=""; ps2pdf  -dEPSCrop ' savePath '.eps '  savePath '.pdf']);
% plot2svg([savePath '.svg'], h);