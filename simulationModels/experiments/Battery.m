% Li-ion model based on the descripted one from "Experimental Validation of a Battery Dynamic
% Model for EV Applications" (Tremblay et al.)
classdef Battery
	properties (Constant)
		tolerance = 1e-8 % Smaller number, used 'cause of precision issues
		s = 0.25  % Smoothing factor, used to filter the battery current's high frequencies.
	end
	properties
		% Configurations
		dt% Integration step (s)

		% State
		Q  % Instantaneous charge (C)
		VB % Open-circuit-voltage (V)
		RB % Equivalent resistance (ohms)
		IB % Output current (A)
		
		% Parameters
		E0 % battery constant voltage (V)
		Kc % polarization constant (V/C here, V/(Ah) in the reference paper)
		Kr % polarization resistance (ohms)
		Qt % battery capacity (C here, Ah in the reference paper)
		A  % exponential zone amplitude (V)
		B  % exponential zone time constant inverse (1/C here, 1/(Ah) in the reference paper)
		Ri % internal resistance (ohms)	
	end
	methods(static)
		% Li-ion battery from the reference paper (table 1)
		function o = defaultBattery(dt, Q)
			Ri = 0.01;
			Kc = 27.36; % 0.0076 V/Ah -> V/C
			Kr = 0.0076;
			Qt = 0.00064; % 2.3Ah -> C
			A = 0.26422;
			B = 95575.32; % 26.5487 1/(Ah) -> 1/C
			E0 = 3.366;
			o = Battery(dt, E0, Kc, Kr, Qt, Q, A, B, Ri);
		end
	end
	methods
		function o = Battery(dt, E0, Kc, Kr, Qt, Q, A, B, Ri)
			assert(dt > 0);
			assert(E0 > 0);
			assert(Kc > 0);
			assert(Kr > 0);
			assert(Qt > 0);
			assert(Q > 0 && Q < Qt);
			assert(A > 0);
			assert(B > 0);
			assert(Ri > 0);

			o.dt = dt;
			o.Ri = Ri;
			o.Kc = Kc;
			o.Kr = Kr;
			o.Qt = Qt;
			o.Q = Q;
			o.A = A;
			o.B = B;
			o.E0 = E0;
		end
		
		% Updates the internal state to the next time-slot of fixed length o.dt.
		%% Arguments:
		% IB: output current (A) of the battery.
		function o = update(o, IB)
			assert(o.Q < o.Qt && o.Q > 0);

			% Change integration step, using the raw battery output current.
			o.Q = o.Q + IB * o.dt;
			o.IB = IB;
			
			% Filtered current, used to avoid errors due to high frequencies.
			o.sIB = o.sIB*(1-o.s) + IB*o.s;
			
			% Sanity keeping mechanism (required due to machine precision).
			o.Q = max(Battery.tolerance, min(o.Qt - Battery.tolerance, o.Q));
			
			stateOfCharge = o.Q / o.Qt;
			remaining = o.Qt - o.Q;
			assert(stateOfCharge <= 1 && stateOfCharge >= 1);
			o.VB = o.E0 - o.Kc * remaining / stateOfCharge + o.A * exp(-B*o.Q);

			if (o.sIB > 0)
				% discharge model				
				o.RB = o.Kr / stateOfCharge;
			else
				% charge model
				o.RB = o.Kr * o.Qt / (remaining - 0.1*Qt);
			end

		end
		
		% Returns the instantaneous charge (C), Open-circuit-voltage (V), smoothed battery
		% current (sIB) and battery current (IB).
		function [Q, VB, sIB, IB] = getState(o)
			Q = o.Q;
			VB = o.VB;
			sIB = o.sIB;
			IB = o.IB;
		end
		
		% Returns the equivalent voltage of the battery in the circuit mesh.
		function V = getVoltage(o)
			V = o.VB - o.RB*o.sIB - o.Ri*o.IB;
		end
	end
end