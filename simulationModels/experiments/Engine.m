classdef Engine
	properties
		% Components
		converter
		load
		battery
	end
	methods(Static)
		function unitTest()
			dt = 0.1; % time increment in secods.
			timeHorizonLen = 1000; % simulation duration in seconds.
			
			power = 5; % load power consumption, W.

			converter = defaultConverter(dt);
			load = staticLoad(power, dt, timeHorizonLen);
			battery = defaultBattery(dt);

			engine = Engine(converter, load, battery);

			% Simulate the whole time horizon with initial currents equal to 0.
			[] = engine.simulate(0,0);
		end
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

			% The initial current to the smooth-mechanism fast initialization.
			o.battery = o.battery.setInitialCurrent(IB);

			while (load.getRemainingIterations() > 0)
				o.load = o.load.update(IL);
				o.battery = o.battery.update(IB);
				[OCV, IL, IB] = o.load.getState();
				[RB, VB] = o.battery.getState();
				o.converter = o.converter.update(IR, RL, VL, RB, VB);
			end
		end
	end
end