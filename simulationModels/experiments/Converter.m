%% Assumed model:
%% -> IR           -> IH           -> (IL-IB)                 <- IB
%%______________INDUCTOR H__________________________________________________
%%        |                      |                    |                     |
%%        |  | (IR-IH)           |  | (IH-IL+IB)      |  | IL               |
%%        |  V                   |  V                 |  V                  |
%%        |                      |                    |                     |
%%   CAPACITOR C1           CAPACITOR C2          RESISTOR RL          RESISTOR RB
%%        |                      |                    |                     |
%%        |                      |                    | +                   | +
%%        |                      |            VOLTAGE SOURCE VL     VOLTAGE SOURCE VB
%%        |                      |                    | -                   | -
%%        |                      |                    |                     |
%%________|______________________|____________________|_____________________|
%%
%% Using Kirchhoff's voltage law:
%% 1: -H(IH-IH0)/dt - (QC2 + dt*(IH-IL+IB))/C2 + (QC1 + dt*(IR-IH))/C1 = 0
%%		<=>  -H/dt*IH + H*IH0/dt - QC2/C2 - dt/C2*IH + dt/C2*IL - dt/C2*IB + QC1/C1 + dt/C1*IR - dt/C1*IH = 0
%%      <=>  (H/dt + dt/C2 + dt/C1)*IH + (-dt/C2)*IL + (dt/C2)*IB = H*IH0/dt - QC2/C2 + QC1/C1 + dt*IR/C1
%% 2: -RL*IL - VL - RB * IB + VB = 0  <=> (0)*IH + (+RL)*IL + (+RB)*IB = VB-VL
%% 3: -RB*IB + VB - (QC2 + dt*(IH-IL+IB))/C2 = 0 <=> -RB*IB + VB - QC2/C2 + dt/C2*IH - dt/C2*IL + dt/C2*IB = 0
%%      <=>  (-dt/C2)*IH + (dt/C2)*IL + (-RB - dt/C2)*IB = QC2/C2 - VB

classdef Converter
	properties(Constant)
		tolerance = 1e-9 % due to machine precision lack
	end
	properties
		% Configurations
		C1 %Capacitance of the first filter capacitor.
		C2 %Capacitance of the second filter capacitor.
		H  %Inductor self-inductance.
		dt %Integration step.
		
		% State
		QC1 %Charge of the first filter capacitor.
		QC2 %Charge of the second filter capacitor.
		IH %Instantaneous inductor current.
		IL %Instantaneous current entering the load.
		IB %Instantaneous current leaving the battery.
	end
	methods(Static)
		function [IR, H, IH, IH0, C1, QC1, C2, QC2, RL, VL, IL, RB, VB, IB, dt] = generateTestInstance()
			dt = rand;
			IR = rand - 0.5;
			H = rand;
			IH = rand - 0.5;
			IH0 = rand - 0.5;
			C1 = rand;
			C2 = rand;
			RL = rand;
			IL = rand;
			RB = rand;
			VB = rand;
			IB = rand;
			
			VL = -RL*IL - RB * IB + VB; % According to (2)
			QC2 = -RB*IB*C2 + VB*C2 - dt*(IH-IL+IB); % According to (3)
			QC1 = H*C1*(IH-IH0)/dt + C1*(QC2 + dt*(IH-IL+IB))/C2 - dt*(IR-IH); % According to (1)
			
			% Sanity checks
			assert(abs(-H*(IH-IH0)/dt - (QC2 + dt*(IH-IL+IB))/C2 + (QC1 + dt*(IR-IH))/C1) < Converter.tolerance);
			assert(abs(-RL*IL - VL - RB * IB + VB) < Converter.tolerance);
			assert(abs(-RB*IB + VB - (QC2 + dt*(IH-IL+IB))/C2) < Converter.tolerance);
		end
		
		function [A, b] = generateModelCoefficients(IR, H, IH0, C1, QC1, C2, QC2, RL, VL, RB, VB, dt)
			A = [(H/dt + dt/C2 + dt/C1), (-dt/C2), (dt/C2);
				 0,                      (+RL),    (+RB);
				 (-dt/C2),               (dt/C2),  (-RB - dt/C2)];
			b = [(H*IH0/dt - QC2/C2 + QC1/C1 + dt*IR/C1);
			     (VB-VL);
				 (QC2/C2 - VB)];
		end
		
		function testModel()
			[IR, H, IH, IH0, C1, QC1, C2, QC2, RL, VL, IL, RB, VB, IB, dt] = Converter.generateTestInstance();
			
			assert(abs((H/dt + dt/C2 + dt/C1)*IH + (-dt/C2)*IL + (dt/C2)*IB -...
				(H*IH0/dt - QC2/C2 + QC1/C1 + dt*IR/C1)) < Converter.tolerance, 'Failed on assertion 1');
			assert(abs(((+RL)*IL + (+RB)*IB) - (VB-VL)) < Converter.tolerance, 'Failed on assertion 2');
			assert(abs((-dt/C2)*IH + (dt/C2)*IL + (-RB - dt/C2)*IB - (QC2/C2 - VB)) < Converter.tolerance,...
				'Failed on assertion 3');
			
			[A, b] = Converter.generateModelCoefficients(IR, H, IH0, C1, QC1, C2, QC2, RL, VL, RB, VB, dt);
			
			assert(max(abs(A*[IH;IL;IB]-b)) < Converter.tolerance, 'Failed on assertion 3');
			
			disp(['SUCCESS for tolerance = ', num2str(Converter.tolerance)]);
		end
	end
	methods
		function o = Converter(C1, C2, H, dt)
			o.C1 = C1;
			o.C2 = C2;
			o.H = H;
			o.dt = dt;
		end											 
		
		% Updates the internal state to the next time-slot of fixed length o.dt.
		%% Arguments:
		% IR: Instantaneous received current.
		% RL: Load current equivalent resistance.
		% VL: Load current equivalent voltage.
		% RB: Battery current equivalent resistance.
		% VB: Battery current equivalent voltage.
		function o = update(o, IR, RL, VL, RB, VB)
			[A, b] = Converter.generateModelCoefficients(IR, H, IH0, C1, QC1, C2, QC2, RL, VL, RB, VB, dt);
			I = A\b;
			o.IH = I(1);
			o.IL = I(2);
			o.IB = I(3);
			o.QC1 = o.QC1 + o.dt * (o.IR - o.IH);
			o.QC2 = o.QC2 + o.dt * (o.IH - o.IL + o.IB);
		end
		
		% Returns the open-circuit-voltage of the converter (V), the current (A) through the load
		% and the output current (A) of the battery.
		function [OCV, IL, IB] = getState(o)
			IL = o.IL;
			IB = o.IB;
			OCV = o.QC1 / o.C1;
		end
	end
end