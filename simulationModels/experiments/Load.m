classdef Load
	properties
		% Configurations
		powerTimeLine
		
		% State
		IL % Input current (A)
		RL % Load resistance (ohms)
		VL % Open-circuit-voltage (V)
		iteration
	end
	methods
		function o = Load()
			o.iteration = 1;
		end
		
		function ri = getRemainingIterations(o)
			length(o.powerTimeLine) - o.iteration + 1;
		end
		
		%% Arguments:
		% IL: Load input current (A)
		function o = update(o, IL)
			assert(o.getRemainingIterations() > 0);
			
			power = o.powerTimeLine(o.iteration);
			
			% updating the state.
			o.IL = IL;
			o.iteration = o.iteration + 1;
		end
		
		% Returns the equivalent load resistance (ohms) and the equivalent "open-circuit-voltage"
		% of the load (V). The actual open-circuit-voltage is zero, since it is a passive device.
		% This is the voltage of the inverse-polarized voltage source used to model the energy
		% consumption together with the load resistor.
		function [RL, VL] = getState(o)
			RL = o.RL;
			VL = o.VL;
		end
	end
end