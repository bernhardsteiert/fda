function axisEqual(posFig)

    aspRatioFig = posFig(3)/posFig(4);
    posSubplot = get(gca,'Position');
    aspRatioSubplot = posSubplot(3)/posSubplot(4);
    set(gca,'YLim',get(gca,'XLim')/(aspRatioFig*aspRatioSubplot))

end