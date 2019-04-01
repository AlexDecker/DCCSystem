%calculate the optimization voltages
function voltage_progression = calc_optimum_voltages(WPTManager,timeSkip)
	voltage_progression = zeros(WPTManager.nt_groups,0);
	for i=1:20
		volt = [];
		for j=1:WPTManager.nt_groups
			volt = [volt;5+10*rand];
		end
		voltage_progression = [voltage_progression,volt];
	end	
end
