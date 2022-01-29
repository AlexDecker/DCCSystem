classdef Battery
	properties (Access = private)
		component
		acquisition
		SOC
	end
	methods
		function bat = Battery(systemName, hierarchy, index)
			% circuit area division
			hierarchy = horizontalCut(hierarchy, [0.5,0.2,0.3]);
			hierarchy.children{3} = verticalCut(hierarchy.children{3}, [0.35,0.3,0.35]);
			
			% building sub-components
			bat.component.name = [systemName, '/bat_', num2str(index)];
			add_block('electricdrivelib/Extra Sources/Battery',...
				bat.component.name,...
				'Position', hierarchy.children{1}.rect ...
			);
			bat.component.hnd = get_param(bat.component.name,'PortHandles');
			set_param(bat.component.name, 'Orientation', 'right');
			
			bat.acquisition.name = [systemName, '/bat_acq_', num2str(index)];
			bat.acquisition.variable = ['bat_acq_', num2str(index)];
			add_block('simulink/Sinks/To Workspace',...
				bat.acquisition.name,...
				'Position', hierarchy.children{3}.children{2}.rect ...
			);
			bat.acquisition.hnd = get_param(bat.acquisition.name,'PortHandles');
			set_param(bat.acquisition.name, 'VariableName', bat.acquisition.variable);
			
			% connecting sub-components
			add_line(systemName,...
				bat.component.hnd.Outport,...
				bat.acquisition.hnd.Inport,...
				'Autorouting', 'on'...
			);
		end
		
		function hnd = positiveHandler(bat)
			hnd = bat.component.hnd.LConn(1);
		end
		
		function hnd = negativeHandler(bat)
			hnd = bat.component.hnd.LConn(2);
		end
		
		function bat = setSOC(bat, SOC)
			set_param(bat.component.name, 'SOC', num2str(SOC));
			bat.SOC = SOC;
		end
		
		function SOC = getSOC(bat)
			SOC = bat.SOC;
		end
	end
end