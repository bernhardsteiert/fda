function setFigure(fh,xfac,yfac,fontsize)

    scaling = 1/(2*fontsize);

    posFig = get(fh,'Position');
    posFig(3) = posFig(3)*xfac;
    posFig(4) = posFig(4)*yfac;
    set(fh,'Position',posFig)
    set(fh,'PaperPosition', [0 0 posFig(3) posFig(4)].*scaling);

end