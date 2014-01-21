function newsite = subplotpos(isite)
    
    if ~mod(ceil(isite/10),2)
        newsite = ceil(isite/10)*10 - mod(isite-1,10);
    else
        newsite = isite;
    end
    
end