function s = siteprop_drug(site)

inh_name = {' + AKTi ',' + MEKi '}; % Cols
ligand_dose = [100 50 20 10 0];
inh_dose = [10 5 2.5 1 0.5 0.25 0]; % Rows

inh_index = ceil(site/10);
dose_index = mod(site-1,10)+1;
if ~mod(inh_index,2)
    dose_index = 11-dose_index;
end
name_ind = ceil(dose_index/5);
dose_index = mod(dose_index-1,5)+1;

s.lig_dose = ligand_dose(dose_index);
s.inh_name = inh_name{name_ind};
s.inh_dose = inh_dose(inh_index);