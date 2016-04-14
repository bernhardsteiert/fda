function d_matrix = plot_inhib_effect()
    cell_lines = {'184A1','MCF10A','SKBR3','HCC1806','HS578T', ...
				'MDA231','BT20','MCF7','T47D'};
	for c = 1:length(cell_lines)
        data = read_data_table(cell_lines{c});
        ligands = unique(data.ligand);
		
        % ERK upon AKT inhibition
		akt_inhib = data((data.meki==0) & (data.akti==1),:);
		akt_noinhib = data((data.meki==0) & (data.akti==0),:);
		
		akt_inhib.diff = akt_inhib.erk_mu - akt_noinhib.erk_mu;
        
		akt_inhib_mean = varfun(@mean, akt_inhib, ...
			'InputVariables', 'diff', ...
			'GroupingVariables', 'ligand');
			
		% AKT upon MEK inhibition
		erk_inhib = data((data.akti==0) & (data.meki==1),:);
		erk_noinhib = data((data.akti==0) & (data.meki==0),:);
		
		erk_inhib.diff = erk_inhib.akt_mu - erk_noinhib.akt_mu;
        
		erk_inhib_mean = varfun(@mean,erk_inhib,...
			'InputVariables', 'diff',...
			'GroupingVariables', 'ligand');
		
		% FOXO upon AKT inhibition
		foxo_akt_inhib = data((data.meki==0) & (data.akti==1),:);
		foxo_akt_noinhib = data((data.meki==0) & (data.akti==0),:);
		
		foxo_akt_inhib.diff = foxo_akt_inhib.foxo_mu - ...
								foxo_akt_noinhib.foxo_mu;
        
		foxo_akt_inhib_mean = varfun(@mean, foxo_akt_inhib, ...
			'InputVariables', 'diff', ...
			'GroupingVariables', 'ligand');
		
		% FOXO upon ERK inhibition
		foxo_erk_inhib = data((data.meki==1) & (data.akti==0),:);
		foxo_erk_noinhib = data((data.meki==0) & (data.akti==0),:);
		
		foxo_erk_inhib.diff = foxo_erk_inhib.foxo_mu - ...
								foxo_erk_noinhib.foxo_mu;
        
		foxo_erk_inhib_mean = varfun(@mean, foxo_erk_inhib, ...
			'InputVariables', 'diff', ...
			'GroupingVariables', 'ligand');
		
		for i = 1:length(ligands)
			dE(c,i) = akt_inhib_mean{ligands{i}, {'mean_diff'}};
			dA(c,i) = erk_inhib_mean{ligands{i}, {'mean_diff'}};
			dFA(c,i) = foxo_akt_inhib_mean{ligands{i}, {'mean_diff'}};
			dFE(c,i) = foxo_erk_inhib_mean{ligands{i}, {'mean_diff'}};
		end
	end
    
	clim = 1;
	
    figure;
    colormap(get_colormap());
    subplot(2,1,1);
    imagesc(dA, [-clim, clim]);
    title('MEKi-driven AKT change')
    colorbar
    set(gca, 'yticklabel', cell_lines);
    set(gca, 'xticklabel', ligands); 
    subplot(2,1,2);
    imagesc(dE, [-clim, clim]);
    title('AKTi-driven ERK change')
    colorbar;
    set(gca, 'yticklabel', cell_lines);
    set(gca, 'xticklabel', ligands);

    
	figure;
    hold on;
	
    subplot(2,1,1);
    title('MEKi driven AKT change')
    hold on;
    plot([0,length(cell_lines)+1],[0,0],'k-');
    plot(dA,'r+');
    xlim([0 length(cell_lines)+1]);
	ylim([-clim clim]);
    set(gca,'xtick',1:length(cell_lines));
    set(gca,'xticklabel',cell_lines)
    set(gca,'xticklabelrotation', 45)

    
    subplot(2,1,2);
    title('AKTi driven ERK change')
    hold on;
    plot([0,length(cell_lines)+1],[0,0],'k-');
    plot(dE,'k+');
    xlim([0 length(cell_lines)+1]);
	ylim([-clim clim]);
    set(gca,'xtick',1:length(cell_lines));
    set(gca,'xticklabel',cell_lines)
    set(gca,'xticklabelrotation', 45)
	
	
	dE_mean = mean(dE,2);
	dA_mean = mean(dA,2);
	dFA_mean = mean(dFA,2);
	dFE_mean = mean(dFE,2);
	d_matrix = [dE_mean, dA_mean, dFA_mean, dFE_mean];
	figure;
	colormap(get_colormap());
	imagesc(d_matrix', [-clim, clim]);
	set(gca, 'xticklabel', cell_lines);
	set(gca,'xticklabelrotation', 90)
	set(gca,'ytick', [1,2,3,4]);
    set(gca, 'yticklabel', {'AKTi->ERK', 'MEKi->AKT', ...
							'AKTi->FOXO3a', 'MEKi->FOXO3a'});
	colorbar;
end

function cm = get_colormap()
	g = [linspace(0.3,1,32)';flipud(linspace(0.3,1,32)')];
	r = [linspace(0.3,1,32)'; ones(32,1)];
    b = [ones(32,1); flipud(linspace(0.3,1,32)')]; 
    cm = [r, g, b];
end
