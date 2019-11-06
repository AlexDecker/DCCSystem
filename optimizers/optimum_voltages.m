%calculate the optimization voltages
function [t,voltage_progression] = calc_optimum_voltages(WPTManager,Pmax,timeSkip)
	voltage_progression = zeros(WPTManager.nt_groups,0);
	for i=1:300000
		volt = [];
		for j=1:WPTManager.nt_groups
			volt = [volt;rand];
		end
		voltage_progression = [voltage_progression,volt];
	end
	[t,voltage_progression] = evalSolution(WPTManager,timeSkip,Pmax,voltage_progression);
	%only the applied voltages remain
	voltage_progression = voltage_progression(:,1:t);
end
