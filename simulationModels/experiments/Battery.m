% Li-ion model based on the descripted one from "Experimental Validation of a Battery Dynamic
% Model for EV Applications" (Tremblay et al.)
% https://www.mdpi.com/2076-3417/8/5/659
classdef Battery
	properties (Constant)
		tolerance = 1e-8 % Smaller number, used 'cause of precision issues
		timeConstant = 30  % secods, used to filter the battery current's high frequencies.
		slidingWindowMinSize = 5 % for the sliding window used in the current filter.
	end
	properties
		% Configurations
		dt % Integration step (h)
		Q0 % Initial charge (Ah)

		% State
		timeSlot % How many states were calculated until now
		extractedCharge % (Ah)
		VB % Equivalent voltage (V)
		RB % Equivalent resistance (ohms) - Unused for now
		outputCurrent % (A)
		filteredOutputCurrent % (A)
		currentSlidingWindow % Sliding window of the outputCurrents used for filtering high frequencies
		
		% Parameters
		E0 % battery voltage constant (V)
		Kc % polarization constant (V/(Ah))
		Kr % polarization resistance (ohms)
		Qt % battery capacity (Ah)
		A  % exponential zone amplitude (V)
		B  % exponential zone time constant inverse (1/(Ah))
		Ri % internal resistance (ohms)	
	end
	methods(Static)
		
		% Seconds to hours conversion.
		function hours = sToHours(s)
			hours = s / 3600;
		end
		
		% Minimal SOC for which this model is expected to work (From 10% to 20% it is expected
		% to work with 10% error for dynamic scenarios).
		function m = minSOC()
			m = Battery.tolerance;
		end
		
		function m = maxSOC()
			m = 100 - Battery.tolerance;
		end
		
		function q = socToCharge(soc, capacity)
			q = capacity * soc / 100;
		end
		
		function soc = chargeToSOC(charge, capacity)
			soc = charge / capacity * 100;
		end

		% Li-ion battery from the reference paper (table 1)
		% dt_seconds: time increment (s).
		% SOC: initial state-of-charge (10-100%).
		% IB: initial output current.
		function o = defaultBattery(dt_seconds, SOC, IB)
			Ri = 0.01;
			Kc = 0.0076;
			Kr = 0.0076;
			Qt = 2.3;
			A = 0.26422;
			B = 26.5487;
			E0 = 3.366;
			o = Battery(dt_seconds, E0, Kc, Kr, Qt, SOC, A, B, Ri, IB);
		end
		
		% time_horizon = Simulation total time in hours.
		% IB = Discharge current, A. Use a negative for charging simulation.
		function simulateForConstantCurrent(battery, time_horizon, IB, showChart)
			dt_hours = battery.dt;
			
			% number of parameter samplings spaced by a (nsamples-1)*ic*dt time horizon.
			nsamples = 46;
			
			maxError = 0.1; % maximum multiplicative variation from the constant current.
			
			% iteration count between samples.
			ic = ceil(time_horizon / (dt_hours * (nsamples-1)));

			Qseries = zeros(nsamples, 1);
			VBseries = zeros(nsamples, 1);
			sIBseries = zeros(nsamples, 1);
			IBseries = zeros(nsamples, 1);
			Itseries = zeros(nsamples, 1);
			
			for j = 1:nsamples-1
				disp(['Processing sample #', num2str(j)]);
				eIB = IB * (1 + (rand-0.5) * maxError);
				[Q, It, VB, sIB, ~] = battery.getState();
				Qseries(j) = Q;
				VBseries(j) = VB;
				sIBseries(j) = sIB;
				IBseries(j) = eIB;
				Itseries(j) = It;
				for i = 1:ic
					battery = battery.update(eIB);
				end
			end
			disp(['Processing sample #', num2str(nsamples)]);
			eIB = IB * (1 + (rand-0.5) * maxError);
			[Q, It, VB, sIB, ~] = battery.getState();
			Qseries(nsamples) = Q;
			VBseries(nsamples) = VB;
			sIBseries(nsamples) = sIB;
			IBseries(j) = eIB;
			Itseries(nsamples) = It;
			
			tAxis = linspace(0, time_horizon, nsamples);
			
			if (showChart)
				figure;
				hold on;
				yyaxis left
				plot(tAxis, VBseries);
				ylabel('(VB)');
				ylim([0 1.5*battery.E0]);
				yyaxis right
				plot(tAxis, Qseries);
				plot(tAxis, Itseries);
				plot(tAxis, IBseries);
				plot(tAxis, sIBseries);
				legend('Battery Voltage','Battery Charge', 'Output Current Integral','Output Current', 'Filtered output current');
				xlabel('Time (h)')
				ylabel('(Ah|A)')
				title('Battery Chart')
				ylim([1.5*min(-abs(IB), battery.minCharge()) 1.5*max(abs(IB), battery.maxCharge())]);
			end
		end
		
		function unitTest(showChart)
			dt_seconds = 0.018; %integration step
			
			%Discharging the battery:
			Q = Battery.maxSOC(); %initial SoC percentage
			IB = 1; %constant discharge current of 1A
			time_horizon = 2.3625; %simulation time in hours
			battery = Battery.defaultBattery(dt_seconds, Q, IB);
			Battery.simulateForConstantCurrent(battery, time_horizon, IB, showChart);
			
			%Charging the battery:
			Q = Battery.minSOC(); %initial SoC percentage
			IB = -1; %constant charge current of 1A
			time_horizon = 2.3625; %simulation time in hours
			battery = Battery.defaultBattery(dt_seconds, Q, IB);
			Battery.simulateForConstantCurrent(battery, time_horizon, IB, showChart);
		end
	end
	methods
	
		function q = minCharge(o)
			q = Battery.socToCharge(Battery.minSOC(), o.Qt);
		end
		
		function q = maxCharge(o)
			q = Battery.socToCharge(Battery.maxSOC(), o.Qt);
		end
		
		function q = instantaneousCharge(o)
			q = o.Qt - o.extractedCharge;
		end
		
		function sanityCheck(o)
			SOC = Battery.chargeToSOC(o.instantaneousCharge(), o.Qt);			
			assert(SOC <= 100 && SOC >= 0);
		end
		
		% Used as sanity keeping mechanism to avoid positive/negative overchage.
		% NOTE: this is a const function. Using the returning values to update the internal state in the
		%	    external scope to avoid rewritting the whole Battery object.
		function extractedCharge = limitExtractedCharge(o)
			minExtractedCharge = o.Qt - o.maxCharge();
			maxExtractedCharge = o.Qt - o.minCharge();
			if (o.extractedCharge > maxExtractedCharge)
				extractedCharge = maxExtractedCharge;
			else
				if (o.extractedCharge < minExtractedCharge)
					extractedCharge = minExtractedCharge;
				else
					extractedCharge = o.extractedCharge;
				end
			end
		end
		
		% Uses a sliding window of the last measured currents. The time correspondence of the
		% window given the integration step is used to control the smoothness of the filteredOutputCurrent
		% result. See the battery datasheet to get the right time constant for the low-pass filter.
		% NOTE: this is a const function. Using the returning values to update the internal state in the
		%	 	external scope to avoid rewritting the whole Battery object.
		function [filteredOutputCurrent, currentSlidingWindow] = updateFilteredCurrent(o)
			oldestCurrent = o.currentSlidingWindow(1);
			windowLength = length(o.currentSlidingWindow);
			currentSlidingWindow = [o.currentSlidingWindow(2:windowLength); o.outputCurrent];
			% Faster mean update, this case is the most repeated one.
			filteredOutputCurrent = o.filteredOutputCurrent + ...
								   (o.outputCurrent - oldestCurrent)/windowLength;
		end
		
		% Creates a new battery model and iterates the first step.
		% dt_seconds: integration step, in seconds.
		% SOC: initial state-of-charge (0-100%)
		% IB: initial output current.
		% The remaining args are the same as specified in attribute block.
		function o = Battery(dt_seconds, E0, Kc, Kr, Qt, SOC, A, B, Ri, IB)
			assert(dt_seconds > 0);
			assert(E0 > 0);
			assert(Kc > 0);
			assert(Kr > 0);
			assert(Qt > 0);
			assert(SOC >= Battery.minSOC() && SOC <= Battery.maxSOC());
			assert(A > 0);
			assert(B > 0);
			assert(Ri > 0);
			
			o.timeSlot = 0;
			o.dt = Battery.sToHours(dt_seconds);
			o.Ri = Ri;
			o.Kc = Kc;
			o.Kr = Kr;
			o.Qt = Qt;
			o.Q0 = Battery.socToCharge(SOC, Qt);
			o.A = A;
			o.B = B;
			o.E0 = E0;
			
			o.extractedCharge = Qt - o.Q0;
			o.extractedCharge = o.limitExtractedCharge();
			
			% Starting as 0 just like the reference paper.
			o.filteredOutputCurrent = 0;

			assert(Battery.slidingWindowMinSize >= 2);
			slidingWindowSize = min(Battery.slidingWindowMinSize, ...
								    ceil(Battery.timeConstant / dt_seconds));
			o.currentSlidingWindow = zeros(slidingWindowSize,1);

			o = o.update(IB);
		end
		
		% Updates the internal state to the next time-slot of fixed length o.dt.
		%% Arguments:
		% IB: output current (A) of the battery.
		function o = update(o, outputCurrent)
			o.sanityCheck();
			
			o.timeSlot = o.timeSlot + 1;
			o.outputCurrent = outputCurrent;

			% Charge integration step, using the raw battery output current.
			o.extractedCharge = o.extractedCharge + o.outputCurrent * o.dt;
			o.extractedCharge = o.limitExtractedCharge();
			o.sanityCheck(); % proc
			
			% Filtered current, used to avoid errors due to high frequencies.
			[o.filteredOutputCurrent, o.currentSlidingWindow] = o.updateFilteredCurrent();
			
			o.VB = o.E0 - o.Ri*o.outputCurrent ...
				 - o.Kc * o.Qt / (o.Qt - o.extractedCharge) * o.extractedCharge ...
				 + o.A * exp(-o.B*o.extractedCharge);
				
			if (o.filteredOutputCurrent > 0)
				% discharge model
				o.VB = o.VB ...
					 - o.Kr * o.Qt / (o.Qt - o.extractedCharge) * o.filteredOutputCurrent;
			else
				% charge model
				o.VB = o.VB ...
					 - o.Kr * o.Qt / (o.extractedCharge - 0.1*o.Qt) * o.filteredOutputCurrent;
			end

		end
		
		% Returns the instantaneous charge (Ah), the extracted charge (Ah),
		% Equivalent voltage (VB), smoothed battery current (sIB) and battery current (IB).
		function [Q, It, VB, sIB, IB] = getState(o)
			o.sanityCheck();
			assert(~isempty(o.extractedCharge));
			assert(~isempty(o.VB));
			assert(~isempty(o.filteredOutputCurrent));
			assert(~isempty(o.outputCurrent));
			Q = o.instantaneousCharge();
			It = o.extractedCharge;
			VB = o.VB;
			sIB = o.filteredOutputCurrent;
			IB = o.outputCurrent;
		end
	end
end