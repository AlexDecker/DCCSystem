classdef Source
	properties (Access = private)
		component
		ammeter
		acquisition
		angularFrequency
		voltage
	end
	methods
		function source = Source(systemName, hierarchy, index)
			% circuit area division
			hierarchy = horizontalCut(hierarchy, [0.25,0.125,0.25,0.125,0.25]);
			hierarchy.children{1} = verticalCut(hierarchy.children{1}, [0.7,0.3]);
			hierarchy.children{3} = verticalCut(hierarchy.children{3}, [0.1,0.25,0.65]);
			hierarchy.children{5} = verticalCut(hierarchy.children{5}, [0.25,0.75]);
			
			% building sub-components
			source.component.name = [systemName, '/src_', num2str(index)];
			add_block('powerlib/Electrical Sources/AC Voltage Source',...
				source.component.name,...
				'Position', hierarchy.children{1}.children{2}.rect ...
			);
			source.component.hnd = get_param(source.component.name,'PortHandles');
			set_param(source.component.name, 'Phase', '0');
			
			source.ammeter.name = [systemName, '/tx_curr_', num2str(index)];
			add_block('powerlib/Measurements/Current Measurement',...
				source.ammeter.name,...
				'Position', hierarchy.children{3}.children{2}.rect ...
			);
			source.ammeter.hnd = get_param(source.ammeter.name,'PortHandles');
			
			source.acquisition.name = [systemName, '/tx_curr_acq_', num2str(index)];
			source.acquisition.variable = ['tx_curr_acq_', num2str(index)];
			add_block('simulink/Sinks/To Workspace',...
				source.acquisition.name,...
				'Position', hierarchy.children{5}.children{1}.rect ...
			);
			source.acquisition.hnd = get_param(source.acquisition.name,'PortHandles');
			set_param(source.acquisition.name, 'VariableName', source.acquisition.variable);
			
			% connecting sub-components
			add_line(systemName,...
				source.component.hnd.RConn,...
				source.ammeter.hnd.LConn,...
				'Autorouting', 'on'...
			);
			
			add_line(systemName,...
				source.ammeter.hnd.Outport,...
				source.acquisition.hnd.Inport,...
				'Autorouting', 'on'...
			);
		end
		
		function hnd = positiveHandler(source)
			hnd = source.ammeter.hnd.RConn;
		end
		
		function hnd = negativeHandler(source)
			hnd = source.component.hnd.LConn;
		end
		
		function source = setVoltage(source, voltage)
			set_param(source.component.name, 'Amplitude', num2str(voltage));
			source.voltage = voltage;
		end
		
		function voltage = getVoltage(source)
			voltage = source.voltage;
		end
		
		function source = setAngularFrequency(source, angularFrequency)
			set_param(source.component.name, 'Frequency', num2str(angularFrequency));
			source.angularFrequency = angularFrequency;
		end
		
		function angularFrequency = getAngularFrequency(source)
			angularFrequency = source.angularFrequency;
		end
	end
end