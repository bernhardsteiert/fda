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
        
    case '01-22-2014'
        celltype = {'MCF10A','MCF10A','MCF10A','MCF10A','MCF10A','MCF10A','184A1','184A1','184A1','184A1','184A1','184A1'};
        ligand_name = {'NS','EPR ','HGF ','EGF ','IGF ','BTC ','BTC ','IGF ','EGF ','HGF ','EPR ','NS'};
        ligand_dose = [0 100 100 100 100 100 100 100 100 100 100 0];
        drug_name = {'Boretezo','Boretezo','Boretezo','Boretezo','No Drug'};
        drug_timing = {'6HR before','1HR before','30MIN after','60MIN after',''};
        
        row = ceil(site/12);
        col = mod(site-1,12)+1;
        if ~mod(row,2)
            col = 13-col;
        end
        
        s.lig_name = ligand_name{col};
        % s.lig_index = lig_index;
        s.lig_dose = ligand_dose(col);
        s.drug_name = drug_name{row};
        s.drug_timing = drug_timing{row};
        s.celltype = celltype{col};
        
    case '01-26-2014'
        ligand_name = {'EGF ','BTC ','EPR ','IGF '};
        ligand_dose = 100;
        drug_name = {' + MEKi ',' + AKTi '};
        drug_dose = [0 0.01 0.1 1];
        celltype = '184A1';
        
        row = ceil(site/12);
        col = mod(site-1,12)+1;
        if ~mod(row,2)
            col = 13-col;
        end
        
        lig_index = ceil(col / 4);
        akti_ind = 5-(mod(col-1,4)+1);
        meki_ind = 5-row;
        if row == 5
            meki_ind = 4 - lig_index + (lig_index < 3);
            lig_index = 4;
        end
        
        s.lig_name = ligand_name{lig_index};
        s.lig_index = lig_index;
        s.lig_dose = ligand_dose;
        s.drug1_name = drug_name{1};
        s.drug1_dose = drug_dose(meki_ind);
        s.drug2_name = drug_name{2};
        s.drug2_dose = drug_dose(akti_ind);
        s.celltype = celltype;

    case '01-27-2014'
        ligand_name = {'EGF ','BTC ','EPR ','IGF '};
        ligand_dose = 100;
        drug_name = {' + MEKi ',' + AKTi '};
        drug_dose = [0 0.01 0.1 1];
        celltype = 'MCF10A';
        
        row = ceil(site/12);
        col = mod(site-1,12)+1;
        if ~mod(row,2)
            col = 13-col;
        end
        
        lig_index = ceil(col / 4);
        akti_ind = 5-(mod(col-1,4)+1);
        meki_ind = 5-row;
        if row == 5
            meki_ind = 4 - lig_index + (lig_index < 3);
            lig_index = 4;
        end
        
        s.lig_name = ligand_name{lig_index};
        s.lig_index = lig_index;
        s.lig_dose = ligand_dose;
        s.drug1_name = drug_name{1};
        s.drug1_dose = drug_dose(meki_ind);
        s.drug2_name = drug_name{2};
        s.drug2_dose = drug_dose(akti_ind);
        s.celltype = celltype;
        
    case '02-02-2014'
        ligand_name = {'EGF ','EGF ','EGF ','EGF ','BTC ','BTC ','BTC ','BTC ','HRG ','HRG ','HRG ','HRG ','NS'};
        lig_index = [1 1 1 1 2 2 2 2 3 3 3 3 4];
        ligand_dose = [100 50 20 10 100 50 20 10 100 50 20 10 0];
        drug_name = ' + MEKi ';
        drug_dose = [10 2.5 2.5/4 2.5/16 0 0];
        celltype = 'MCF10A';
        
        row = ceil(site/12);
        col = mod(site-1,12)+1;
        if ~mod(row,2)
            col = 13-col;
        end
        
        if row == 6
            col = 13;
        end
        
        s.lig_name = ligand_name{col};
        s.lig_index = lig_index(col);
        s.lig_dose = ligand_dose(col);
        s.drug_name = drug_name;
        s.drug_dose = drug_dose(row);
        s.celltype = celltype;
        
    case '02-12-2014-wtAkt'
        ligand_name = {'BTC','EPR','EGF','HGF','HRG','IGF','NS'};
        celltype = {'Native', 'ERKmut'};
        ligand_dose = [6.25 12.5 25 50 100 100 50 25 12.5 6.25 0];
        
        row = ceil(site/10);
        col = mod(site-1,10)+1;
        if ~mod(row,2)
            col = 11-col;
        end
        if row == 7
            col = 11;
        end
        
        s.lig_name = ligand_name{row};
        s.lig_index = row;
        s.lig_dose = ligand_dose(col);
        s.celltype = celltype{(col>5)+1};
        
    case '02-26-2014'
        ligand_name = {'BTC','EPR','EGF','HGF','HRG','IGF','NS'};
        celltype = {'Native', 'ERKmut'};
        ligand_dose = [6.25 12.5 25 50 100 100 50 25 12.5 6.25 0];
        
        row = ceil(site/10);
        col = mod(site-1,10)+1;
        if ~mod(row,2)
            col = 11-col;
        end
        if row == 7
            col = 11;
            celltype{2} = celltype{1};
        end
        
        s.lig_name = ligand_name{row};
        s.lig_index = row;
        s.lig_dose = ligand_dose(col);
        s.celltype = celltype{(col>5)+1};
        
    case {'02-15-2014','02-15-2014_retracked','02-15-2014_retracked_nucl'}
        ligand_name = {'BTC','EPR','EGF','HGF','HRG','IGF','NS'};
        ligand_dose = [4 20 100 4 20 100 4 20 100 0];
        sensor_name = {'EKAREV','EKAREV','EKAREV','FOXO3a','FOXO3a','FOXO3a','Dual_EKAREV_FOXO3a','Dual_EKAREV_FOXO3a','Dual_EKAREV_FOXO3a'};
        
        row = ceil(site/9);
        col = mod(site-1,9)+1;
        if ~mod(row,2)
            col = 10-col;
        end
        col_bak = col;
        if row == 7
            col = 10;
        end
        
        s.lig_name = ligand_name{row};
        s.lig_index = row;
        s.lig_dose = ligand_dose(col);
        s.sensor_name = sensor_name{col_bak};
        
        
    case '02-27-2014'
        ligand_name = {'BTC','EPR','EGF','HGF','HRG','IGF','NS'};
        ligand_dose = [4 20 100 4 20 100 4 20 100 0];
        sensor_name = {'EKAREV','EKAREV','EKAREV','FOXO3a','FOXO3a','FOXO3a','Dual_EKAREV_FOXO3a','Dual_EKAREV_FOXO3a','Dual_EKAREV_FOXO3a'};
        
        row = ceil(site/9);
        col = mod(site-1,9)+1;
        if ~mod(row,2)
            col = 10-col;
        end
        col_bak = col;
        if row == 7
            col = 10;
        end
        
        s.lig_name = ligand_name{row};
        s.lig_index = row;
        s.lig_dose = ligand_dose(col);
        s.sensor_name = sensor_name{col_bak};
        
    case {'03-30-2014','03-30-2014_cleaned'}
        igf_name = 'IGF';
        igf_dose = [.8 4 20 100 .8 4 20 100 .8 4 20 100];
        ligand_name = {'NS','BTC','EGF','EGF','NS','BTC','EGF','EGF','NS','BTC','EGF','EGF'};
        ligand_dose = [0 100 100 1 0 100 100 1 0 100 100 1];
        celltype = {'MCF10A','MCF10A','MCF10A','MCF10A','184A1','184A1','184A1','184A1','HCC1806','HCC1806','HCC1806','HCC1806'};
        drug_name = 'AKTi';
        drug_dose = [.5 .1 .025 .00625 0 0];
        
        row = ceil(site/12);
        col = mod(site-1,12)+1;
        if ~mod(row,2)
            col = 13-col;
        end
        
        if row == 6
            s.lig_name = ligand_name{col};
            s.lig_dose = ligand_dose(col);
        else
            s.lig_name = igf_name;
            s.lig_dose = igf_dose(col);
        end
        
        s.celltype = celltype{col};
        s.drug_name = drug_name;
        s.drug_dose = drug_dose(row);
        
    case {'04-04-2014','04-04-2014_all_cleaned'}
        ligand_name = {'EGF','EGF','EGF','EGF','IGF','IGF','IGF','IGF'};
        celltype = {'MCF10A','MCF10A','MCF10A','MCF10A','MCF10A','184A1','184A1','184A1','184A1','184A1'};
        meki_dose = [0 1/6^3 1/6^2 1/6 1 0 1/6^3 1/6^2 1/6 1];
        akti_dose = [0 .005 .05 .5 0 .005 .05 .5];
        
        row = ceil(site/10);
        col = mod(site-1,10)+1;
        if ~mod(row,2)
            col = 11-col;
        end
        
        s.lig_name = ligand_name{row};
        s.celltype = celltype{col};
        s.drug1_name = 'MEKi';
        s.drug2_name = 'AKTi';
        s.drug1_dose = meki_dose(col);
        s.drug2_dose = akti_dose(row);
        
    case '03-23-2014-10A1806'
        ligand_name = {'EGF','EGF','EGF','IGF','IGF','IGF'};
        ligand_dose = [.1 .1 .1 1 1 1 1 1 1 10 10 10];
        celltype = {'MCF10A','MCF10A','MCF10A','MCF10A','MCF10A','MCF10A','HCC1806','HCC1806','HCC1806','HCC1806','HCC1806','HCC1806'};
        meki_dose = [0 0 1 1 1 1 1 1 1 1 0 0];
        meki_timing = {'No MEKi','No MEKi','Late (2HR)','Late (2HR)','Early (30MIN)','Early (30MIN)','Early (30MIN)','Early (30MIN)','Late (2HR)','Late (2HR)','No MEKi','No MEKi'};
        akti_dose = [0 .01 .1 0 .01 .1];
        
        row = ceil(site/12);
        col = mod(site-1,12)+1;
        if ~mod(row,2)
            col = 13-col;
        end
        
        s.lig_name = ligand_name{row};
        s.lig_dose = ligand_dose(row+mod(col-1,2)*6);
        s.celltype = celltype{col};
        s.drug1_name = 'MEKi';
        s.drug2_name = 'AKTi';
        s.drug1_dose = meki_dose(col);
        s.drug2_dose = akti_dose(row);
        s.drug1_timing = meki_timing{col};
        
    case '04-15-2014'
        celltype = {'MCF10A-WT','MCF10A-WT','MCF10A-AKTspec','MCF10A-AKTspec','MCF10A-ERKspec','MCF10A-ERKspec','184A1-WT','184A1-WT','184A1-AKTspec','184A1-AKTspec','184A1-ERKspec','184A1-ERKspec'};
        ligand_name = {'IGF','HRG','HGF','EGF','BTC','NS'};
        ligand_dose = [100 20 100 20 100 20 100 20 100 20 100 20];
        
        row = ceil(site/12);
        col = mod(site-1,12)+1;
        if ~mod(row,2)
            col = 13-col;
        end
        
        s.celltype = celltype{col};
        s.lig_name = ligand_name{row};
        s.lig_index = row;
        s.lig_dose = ligand_dose(col);
        if row == 6
            s.celltype = celltype{col-(12-col)};
            s.lig_name = ligand_name{row};
            s.lig_dose = 0;
        end
        
    case '04-18-2014'
        celltype = {'MCF10A','MCF10A','MCF10A','MCF10A','MCF10A','MCF10A','184A1','184A1','184A1','184A1','184A1','184A1'};
        ligand_name = 'EGF';
        ligand_dose = [100 20 4 .8 .16 0];
        drug_name = 'MEKi';
        drug_dose = [.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0 0 .1/4^4 .1/4^3 .1/4^2 .1/4 .1];
        
        row = ceil(site/12);
        col = mod(site-1,12)+1;
        if ~mod(row,2)
            col = 13-col;
        end
        
        s.celltype = celltype{col};
        s.lig_name = ligand_name;
        s.lig_dose = ligand_dose(row);
        s.drug_name = drug_name;
        s.drug_dose = drug_dose(col);
        
    case '03-24-2014'
        celltype = 'MCF10A';
        ligand_name = {'IGF','HRG','HGF','EGF','BTC','EPR'};
        ligand_dose = [100 100/2.5 100/2.5^2 100/2.5^3 100/2.5^4 100/2.5^5 100/2.5^6 100/2.5^7 100/2.5^8 0];
        
        row = ceil(site/10);
        col = mod(site-1,10)+1;
        if ~mod(row,2)
            col = 11-col;
        end
        
        s.celltype = celltype;
        s.lig_name = ligand_name{row};
        s.lig_dose = ligand_dose(col);
        
    case {'02-20-2015','02-20-2015_cleaned'}
        
        meki_dose = [.125/3^0 .125/3^1 .125/3^2 .125/3^3 .125/3^4 0]; % Rows
        akti_dose = [1/1.5^0 1/1.5^1 1/1.5^2 1/1.5^3 1/1.5^4 1/1.5^5 1/1.5^6 1/1.5^7 1/1.5^8 1/1.5^9 1/1.5^10 0]; % Cols

        row = ceil(site/12);
        col = mod(site-1,12)+1;
        if ~mod(row,2)
            col = 13-col;
        end
        
        s.lig_name = 'EGF';
        s.lig_dose = 20;
        s.drug1_name = 'AKTi';
        s.drug2_name = 'MEKi';
        s.drug1_dose = akti_dose(col);
        s.drug2_dose = meki_dose(row);
         
    otherwise
        error('Unknown data-set!!')
        
end