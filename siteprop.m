function s = siteprop(site)

ligand_name = {'EGF ','IGF ','FGF ','HRG ','HGF ','EPR ','BTC '}; % Rows
inh_name = {' + MEKi',' + AKTi',' + PI3Ki','','','','','','',''}; % Cols
ligand_dose = [100 100 100 100 50 20 10 5 2.5 0];
inh_dose = [10 10 10 0 0 0 0 0 0 0];

lig_index = ceil(site/10);
inh_index = mod(site-1,10)+1;
if ~mod(lig_index,2)
    inh_index = 11-inh_index;
end

s.lig_name = ligand_name{lig_index};
s.lig_index = lig_index;
s.lig_dose = ligand_dose(inh_index);
s.inh_name = inh_name{inh_index};
s.inh_dose = inh_dose(inh_index);
s.col_index = inh_index;