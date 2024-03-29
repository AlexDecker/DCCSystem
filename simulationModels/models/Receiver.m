classdef Receiver
	properties (Access = private)
		ammeter
        rc
        rectifier
        load
        battery
	end
	methods
		function obj = Receiver(systemName, hierarchy, index)
			% circuit area division
			hierarchy = horizontalCut(hierarchy, [0.23,0.07,0.13,0.04,0.23,0.04,0.3]);
            hierarchy.children{1} = verticalCut(hierarchy.children{1}, [0.5,0.1,0.4]);
            hierarchy.children{3} = verticalCut(hierarchy.children{3}, [0.2,0.8]);
            hierarchy.children{5} = verticalCut(hierarchy.children{5}, [0.45,0.3,0.25]);
            hierarchy.children{7} = verticalCut(hierarchy.children{7}, [0.2,0.8]);
			
			% building sub-components
			obj.ammeter = RXAmmeter(systemName, hierarchy.children{1}.children{1}, index);
            obj.rc = RCBranch(systemName, hierarchy.children{1}.children{3}, 'rx', index);
            obj.rectifier = Rectifier(systemName, hierarchy.children{3}.children{2}, index);
            obj.load = Load(systemName, hierarchy.children{5}.children{2}, index);
            obj.battery = Battery(systemName, hierarchy.children{7}.children{2}, index);

			% connecting sub-components %%%%%%%%%%%
            
            % The battery is connected using a positive-to-positive and negative-to-negative
            % layout since the input current is meant to recharge it.
			add_line(systemName,...
				obj.rectifier.inNegativeHandler(),...
				obj.battery.negativeHandler(),...
				'Autorouting', 'on'...
			);
            add_line(systemName,...
				obj.rectifier.inPositiveHandler(),...
				obj.battery.positiveHandler(),...
				'Autorouting', 'on'...
			);
            
            add_line(systemName,...
				obj.rectifier.inNegativeHandler(),...
				obj.load.positiveHandler(),...
				'Autorouting', 'on'...
			);
            add_line(systemName,...
				obj.rectifier.inPositiveHandler(),...
				obj.load.negativeHandler(),...
				'Autorouting', 'on'...
			);
            
            add_line(systemName,...
				obj.rectifier.outNegativeHandler(),...
				obj.ammeter.positiveHandler(),...
				'Autorouting', 'on'...
			);
            add_line(systemName,...
				obj.rectifier.outPositiveHandler(),...
				obj.rc.negativeHandler(),...
				'Autorouting', 'on'...
			);
			
		end
		
		function hnd = positiveHandler(obj)
			hnd = obj.rc.positiveHandler();
		end
		
		function hnd = negativeHandler(obj)
			hnd = obj.ammeter.negativeHandler();
		end
		
		function obj = setResistance(obj, resistance)
			obj.rc = obj.rc.setResistance(resistance);
		end
		
		function resistance = getResistance(obj)
			resistance = obj.rc.getResistance();
		end
		
		function obj = setCapacitance(obj, capacitance)
			obj.rc = obj.rc.setCapacitance(capacitance);
		end
		
		function obj = setConsummerResistance(obj, resistance)
			obj.load = obj.load.setResistance(resistance);
		end
		
		function resistance = getConsummerResistance(obj)
			resistance = obj.load.getResistance();
		end
		
		function obj = setSOC(obj, SOC)
			obj.battery = obj.battery.setSOC(SOC);
		end
		
		function SOC = getSOC(obj)
			SOC = obj.battery.getSOC();
		end
		
		% get acquisition object names %
		
		function acqName = getLoadCurrentAcquisition(obj)
			acqName = obj.load.acquisition.variable;
		end
		
		function acqName = getReceivingCurrentAcquisition(obj)
			acqName = obj.ammeter.acquisition.variable;
		end
		
		function acqName = getReceivingVoltageAcquisition(obj)
			acqName = obj.rc.acquisition.variable;
		end
		
		function acqName = getBatteryAcquisition(obj)
			acqName = obj.battery.acquisition.variable;
		end
		
	end
end