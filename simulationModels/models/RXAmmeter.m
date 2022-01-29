classdef RXAmmeter
	properties (Access = private)
		ammeter
		acquisition
	end
	methods
		function obj = RXAmmeter(systemName, hierarchy, index)
			% circuit area division
			hierarchy = horizontalCut(hierarchy, [0.4,0.2,0.4]);
			hierarchy.children{1} = verticalCut(hierarchy.children{1}, [0.6,0.4]);
			hierarchy.children{3} = verticalCut(hierarchy.children{3}, [0.4,0.6]);
			
			% building sub-components
			obj.ammeter.name = [systemName, '/rx_curr_', num2str(index)];
			add_block('powerlib/Measurements/Current Measurement',...
				obj.ammeter.name,...
				'Position', hierarchy.children{3}.children{2}.rect ...
			);
			obj.ammeter.hnd = get_param(obj.ammeter.name,'PortHandles');
			set_param(obj.ammeter.name,'Orientation', 'left');
			
			obj.acquisition.name = [systemName, '/rx_curr_acq_', num2str(index)];
			obj.acquisition.variable = ['rx_curr_acq_', num2str(index)];
			add_block('simulink/Sinks/To Workspace',...
				obj.acquisition.name,...
				'Position', hierarchy.children{1}.children{1}.rect ...
			);
			obj.acquisition.hnd = get_param(obj.acquisition.name,'PortHandles');
			set_param(obj.acquisition.name, 'VariableName', obj.acquisition.variable);
			set_param(obj.acquisition.name,'Orientation', 'left');
			
			% connecting sub-components
			add_line(systemName,...
				obj.ammeter.hnd.Outport,...
				obj.acquisition.hnd.Inport,...
				'Autorouting', 'on'...
			);
			
		end
		
		function hnd = positiveHandler(obj)
			hnd = obj.ammeter.hnd.LConn;
		end
		
		function hnd = negativeHandler(obj)
			hnd = obj.acquisition.hnd.RConn;
		end
	end
end