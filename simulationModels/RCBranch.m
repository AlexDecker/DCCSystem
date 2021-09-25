classdef RCBranch
	properties (Access = private)
		component
		voltmeter
		acquisition
		capacitance
		resistance
	end
	methods
		function [rc, hierarchy] = RCBranch(systemName, hierarchy, prefix, index)
			% circuit area division
			hierarchy = verticalCut(hierarchy, [0.33,0.33,0.34]);
			hierarchy.children{1} = horizontalCut(hierarchy.children{1}, [0.2,0.6,0.2]);
			hierarchy.children{3} = horizontalCut(hierarchy.children{3}, [0.2,0.3,0.2,0.3]);
			
			% building sub-components
			rc.component.name = [systemName, '/', prefix, '_rc_', num2str(index)];
			add_block('powerlib/Elements/Series RLC Branch',...
				rc.component.name,...
				'Position', hierarchy.children{1}.children{2}.rect ...
			);
			rc.component.hnd = get_param(rc.component.name,'PortHandles');
			set_param(rc.component.name, 'BranchType', 'RC');
			
			rc.voltmeter.name = [systemName, '/', prefix, '_volt_', num2str(index)];
			add_block('powerlib/Measurements/Voltage Measurement',...
				rc.voltmeter.name,...
				'Position', hierarchy.children{3}.children{2}.rect ...
			);
			rc.voltmeter.hnd = get_param(rc.voltmeter.name,'PortHandles');
			
			rc.acquisition.name = [systemName, '/', prefix, '_volt_acq_', num2str(index)];
			rc.acquisition.variable = [prefix, '_volt_acq_', num2str(index)];
			add_block('simulink/Sinks/To Workspace',...
				rc.acquisition.name,...
				'Position', hierarchy.children{3}.children{4}.rect ...
			);
			rc.acquisition.hnd = get_param(rc.acquisition.name,'PortHandles');
			set_param(rc.acquisition.name, 'VariableName', rc.acquisition.variable);
			
			% connecting sub-components
			add_line(systemName,...
				rc.component.hnd.LConn,...
				rc.voltmeter.hnd.LConn(1),...
				'Autorouting', 'on'...
			);
			
			add_line(systemName,...
				rc.component.hnd.RConn,...
				rc.voltmeter.hnd.LConn(2),...
				'Autorouting', 'on'...
			);
			
			add_line(systemName,...
				rc.voltmeter.hnd.Outport,...
				rc.acquisition.hnd.Inport,...
				'Autorouting', 'on'...
			);
		end
		
		function hnd = positiveHandler(rc)
			hnd = rc.component.hnd.LConn;
		end
		
		function hnd = negativeHandler(rc)
			hnd = rc.component.hnd.RConn;
		end
		
		function rc = setResistance(rc, resistance)
			set_param(rc.component.name, 'Resistance', num2str(resistance));
			rc.resistance = resistance;
		end
		
		function resistance = getResistance(rc)
			resistance = rc.resistance;
		end
		
		function rc = setCapacitance(rc, capacitance)
			set_param(rc.component.name, 'Capacitance', num2str(capacitance));
			rc.capacitance = capacitance;
		end
		
		function capacitance = getCapacitance(rc)
			capacitance = rc.capacitance;
		end
	end
end