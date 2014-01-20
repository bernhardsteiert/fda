function s = siteprop(site,dataset)

if(~exist('dataset','var'))
	dataset = 'default';
end

switch dataset
    case 'default'

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
        
    case '2D_dose_response_drugsVSEGF_130903'
        
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
        
        s.lig_name = 'EGF';
        s.lig_dose = ligand_dose(dose_index);
        s.inh_name = inh_name{name_ind};
        s.inh_ind = name_ind;
        s.inh_dose = inh_dose(inh_index);
        
    case '2D_dose_response_drugsVSHGF_130826'
        
        error('New data set, check and copy over siteprop from 2D_dose_response_drugsVSEGF_130903 ...')
        
    case '12-08-2013'
        
        ligand_dose = [1 10 100 0]; % Rows
        ligand_name = {'FGF ','HRG ','IGF ','EGF ','EGF ','BTC ','BTC ','HGF ','EPR ','NS'}; % Cols
        celltype = {'Native', 'ERKmut'};
        
        dose_index = ceil(site/10);
        lig_index = mod(site-1,10)+1;
        if ~mod(dose_index,2)
            lig_index = 11-lig_index;
        end
        name_ind = 1+(site>30);
        if dose_index > 3
            dose_index = 7-dose_index;
        end
        if lig_index == 10
            dose_index = 4;
        end
        
        s.lig_name = ligand_name{lig_index};
        s.lig_index = lig_index;
        s.lig_dose = ligand_dose(dose_index);
        s.celltype = celltype{name_ind};
        
    case '12-25-2013'
        
        ligand_name = {'NS','BTC ','EGF ','IGF '}; % Rows (1st / Last Col)
        ligand_dose = [0 1 10 100]; % Cols
        inh_dose = [0 1];
        inh_name = {'',' + MEKi'}; % Cols
        celltype = {'Native', 'ERKmut'};
        
        lig_index = mod(ceil(site/10)-1,3)+2;
        dose_index = mod(site-1,10)+1;
        if mod(ceil(site/10)+1,2)
            dose_index = 11-dose_index;
        end
        inh_index = ceil((dose_index)/5);
        if inh_index == 2
            dose_index = 11-dose_index;
        end
        if dose_index == 1
            lig_index = 1;
        elseif dose_index > 3 && dose_index < 8
            dose_index = 4;
        end
        name_ind = 1+(site>30);
        
        s.lig_name = ligand_name{lig_index};
        s.lig_index = lig_index;
        s.lig_dose = ligand_dose(dose_index);
        s.inh_name = inh_name{inh_index};
        s.inh_dose = inh_dose(inh_index);
        s.celltype = celltype{name_ind};
        
    otherwise
        error('Unknown data-set!!')
        
end