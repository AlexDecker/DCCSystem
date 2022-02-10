classdef Transmitter
	properties (Access = private)
		source
        rc
	end
	methods
		function obj = Transmitter(systemName, hierarchy, index)
			% circuit area division
			hierarchy = verticalCut(hierarchy, [0.5,0.1,0.4]);
            hierarchy.children{1} = horizontalCut(hierarchy.children{1}, [0.8,0.2]);
            hierarchy.children{3} = horizontalCut(hierarchy.children{3}, [0.2,0.8]);
			
			% building sub-components
			obj.source = Source(systemName, hierarchy.children{1}.children{1}, index);
			obj.rc = RCBranch(systemName, hierarchy.children{3}.children{2}, 'tx', index);

			% connecting sub-components
			add_line(systemName,...
				obj.source.negativeHandler(),...
				obj.rc.positiveHandler(),...
				'Autorouting', 'on'...
			);
			
		end
		
		function hnd = positiveHandler(obj)
			hnd = obj.source.positiveHandler();
		end
		
		function hnd = negativeHandler(obj)
			hnd = obj.rc.negativeHandler();
		end
		
		function obj = setVoltage(obj, voltage)
			obj.source = obj.source.setVoltage(voltage);
		end
		
		function voltage = getVoltage(obj)
			voltage = obj.source.getVoltage();
		end
		
		function obj = setFrequency(obj, frequency)
			obj.source = obj.source.setFrequency(frequency);
		end
		
		function frequency = getFrequency(obj)
			frequency = obj.source.getFrequency();
		end
		
		function obj = setResistance(obj, resistance)
			obj.rc = obj.rc.setResistance(resistance);
		end
		
		function obj = setCapacitance(obj, capacitance)
			obj.rc = obj.rc.setCapacitance(capacitance);
		end
		
		function resistance = getResistance(obj)
			resistance = obj.rc.getResistance();
		end
	end
end