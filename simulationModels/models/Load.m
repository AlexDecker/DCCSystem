classdef Load
	properties (SetAccess = private)
		component
		ammeter
		acquisition
		resistance
	end
	methods
		function load = Load(systemName, hierarchy, index)
			% circuit area division
			hierarchy = horizontalCut(hierarchy, [0.15,0.05,0.3,0.2,0.3]);
			hierarchy.children{1} = verticalCut(hierarchy.children{1}, [0.2,0.3,0.5]);
			hierarchy.children{3} = verticalCut(hierarchy.children{3}, [0.6,0.3,0.1]);
			hierarchy.children{5} = verticalCut(hierarchy.children{5}, [0.3,0.3,0.4]);
			
			% building sub-components
			load.component.name = [systemName, '/load_rc_', num2str(index)];
			add_block('powerlib/Elements/Series RLC Branch',...
				load.component.name,...
				'Position', hierarchy.children{1}.children{2}.rect ...
			);
			load.component.hnd = get_param(load.component.name,'PortHandles');
			set_param(load.component.name, 'BranchType', 'R');
			set_param(load.component.name, 'Orientation', 'down');
			
			load.ammeter.name = [systemName, '/load_curr_', num2str(index)];
			add_block('powerlib/Measurements/Current Measurement',...
				load.ammeter.name,...
				'Position', hierarchy.children{3}.children{2}.rect ...
			);
			load.ammeter.hnd = get_param(load.ammeter.name,'PortHandles');
			
			load.acquisition.name = [systemName, '/load_curr_acq_', num2str(index)];
			load.acquisition.variable = ['load_curr_acq_', num2str(index)];
			add_block('simulink/Sinks/To Workspace',...
				load.acquisition.name,...
				'Position', hierarchy.children{5}.children{2}.rect ...
			);
			load.acquisition.hnd = get_param(load.acquisition.name,'PortHandles');
			set_param(load.acquisition.name, 'VariableName', load.acquisition.variable);
			
			% connecting sub-components
			add_line(systemName,...
				load.component.hnd.RConn,...
				load.ammeter.hnd.LConn,...
				'Autorouting', 'on'...
			);
			
			add_line(systemName,...
				load.ammeter.hnd.Outport,...
				load.acquisition.hnd.Inport,...
				'Autorouting', 'on'...
			);
		end
		
		function hnd = positiveHandler(load)
			hnd = load.ammeter.hnd.RConn;
		end
		
		function hnd = negativeHandler(load)
			hnd = load.component.hnd.LConn;
		end
		
		function load = setResistance(load, resistance)
			set_param(load.component.name, 'Resistance', num2str(resistance));
			load.resistance = resistance;
		end
		
		function resistance = getResistance(load)
			resistance = load.resistance;
		end
	end
end