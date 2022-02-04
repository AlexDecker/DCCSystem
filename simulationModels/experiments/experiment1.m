% This script aims on gather data from a 1TX1RX wpt circuit, in order to model
% the set converter/battery/consumer in terms only of the state-of-charge and the
% consuming current (considering the other parameters of the battery are constant)

generateData = true;

if generateData
	wpt = WPTSystem('wpt', 1, 1, [], 'experiments\ParameterGeneration\');

	% Change voltages, resistances, couplings and SOC
	wpt = wpt.updateParameters();

	wpt.run();

	settings = wpt.getSettings();
	
	disp('Input Settings:');
	disp(['Voltage (V): ', num2str(settings.voltages(1)),...
		  ' TX Impedance (ohms): ', num2str(settings.internalImpedance(1)),...
		  ' RX Fixed Impedance (ohms): ', num2str(settings.internalImpedance(2))...
		  );
	disp(['SOC (0-1): ', num2str(settings.SOC(1)),...
		  ' MutualInductance (H): ', num2str(settings.M(1,2)),...
		  ' Frequency (Hz): ', num2str(settings.f),...
		  ' Consumer Resistance (ohms): ', num2str(settings.RL(1)),...
		  );
	
	% To process signal readings and calculate the equivalent impedance
	resMgr = ResistanceAbstractionManager(settings.f);
	
	

	wpt.destroy();
else
	% Analyze data
end