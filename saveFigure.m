function saveFigure(h, name, renameAxis)

if(~exist('renameAxis','var'))
	renameAxis = true;
elseif isempty(renameAxis)
    renameAxis = true;
end

savePath =  ['./Figures/' datestr(now, 30) '_' name];

if renameAxis
    axis_handle = get(gcf, 'CurrentAxes');
    xlabel(axis_handle, 'time [min]')
    ylabel(axis_handle, 'log10(FOXO3a Ratio [C/N])')
end

saveas(h, savePath, 'fig');
print('-depsc', savePath);
% print('-dpng', savePath);
system(['export LD_LIBRARY_PATH=""; ps2pdf  -dEPSCrop ' savePath '.eps '  savePath '.pdf']);
% plot2svg([savePath '.svg'], h);