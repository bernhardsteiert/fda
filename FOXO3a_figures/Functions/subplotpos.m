function newsite = subplotpos(isite,nCol)
    if(~exist('nCol','var'))
        nCol = 10;
    end
    
    if ~mod(ceil(isite/nCol),2)
        newsite = ceil(isite/nCol)*nCol - mod(isite-1,nCol);
    else
        newsite = isite;
    end
    
end