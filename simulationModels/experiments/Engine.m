classdef Engine
	properties
		% Components
		converter
		load
		battery
	end
	methods
		function o = Engine(converter, load, battery)
			o.converter = converter;
			o.load = load;
			o.battery = battery;
		end											 
		
		%% Arguments:
		% IL: Initial load current (A)
		% IB: Initial battery output current (A)
		function [] = simulate(o, IL, IB)
			while (load.getRemainingIterations() > 0)
				o.load = o.load.update(IL);
				o.battery = o.battery.update(IB);
				[RL, VL] = o.load.getState();
				[RB, VB] = o.battery.getState();
				o.converter = o.converter.update(IR, RL, VL, RB, VB);
			end
		end
	end
end