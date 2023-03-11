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
	methods(static)
		function o = staticLoad(power, dt, timeHorizonLength)
			
		end
		
		% Creates a random timeline based on user specifications.
		% minPower: Minimal power within the time horizon, in watts.
		% maxPower: Maximal power within the time horizon, in watts.
		% dt: time increment of the simulation.
		% timeHorizonLength: time duration of the simulation.
		% nPoints: number of random points before interpolation. Minimum of 4. Controls
		%		   the smoothness of the curve.
		function o = randomLoad(minPower, maxPower, dt, timeHorizonLength, nPoints)
			assert(minPower >= 0);
			assert(maxPower >= minPower);
			assert(dt > 0);
			assert(timeHorizonLength > dt);
			assert(nPoints >= 4);
			powerTimeLine = zeros(nPoints, 2);
			powerTimeLine(:, 1) = rand(nPoints, 1
			o = Load(powerTimeLine, dt);
		end
		
		function o = defaultLoad(dt)
		end
	end
	methods
		% powerTimeLine: Lookup table where the first col is the time and the second is
		%				 the dissipated power of the load in that time. Must have at least
		%				 rows.
		% dt: time increment of the simulation.
		function o = Load(powerTimeLine, dt)
			o.iteration = 1;
			o.RL = 0; % for now, we are just considering the power source abstraction.
			assert(dt > 0);
			[rows, cols] = size(powerTimeLine);
			assert(rows >= 4);
			assert(cols == 2);
			% sort by time
			powerTimeLine = sortrows(powerTimeLine, 1);
			% homogeneous domain
			tDomain = powerTimeLine(1,1):dt:powerTimeLine(rows,1);
			% finding the powerTimeLine vector using shape-preserving cubic interpolation.
			o.powerTimeLine = interp1(powerTimeLine(:,1),powerTimeLine(:,2),tDomain,'pchip')
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