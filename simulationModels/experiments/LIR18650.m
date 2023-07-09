% LIR18650 Li-Ion battery data to compare with the results of our battery model.
classdef LIR18650
	properties (Constant)
		dischargeVoltage1C = [4.0795, 3.9099, 3.8040, 3.6556, 3.6026, 3.5179, 3.4649, 3.4225, 3.3589, 3.2742, 3.1894, 3.0199, 2.8079]
		dischargeSOC1C = [0, 0.0864, 0.2006, 0.4012, 0.5000, 0.7377, 0.8519, 0.8951, 0.9259, 0.9475, 0.9630, 0.9784, 0.9907]
		dischargeVoltage0_5C = [4.1430, 4.0053, 3.8887, 3.7510, 3.7086, 3.6238, 3.5921, 3.5603, 3.5073, 3.4543, 3.3695, 3.2318, 3.1470, 3.0411]
		dischargeSOC0_5C = [0, 0.0864, 0.2006, 0.4012, 0.5000, 0.7377, 0.8519, 0.8951, 0.9259, 0.9475, 0.9630, 0.9784, 0.9907, 1.0000]
		dischargeVoltage0_2C = [4.1642, 4.0477, 3.9417, 3.8146, 3.7722, 3.6980, 3.6450, 3.6238, 3.6132, 3.5603, 3.5073, 3.3801, 3.2848, 3.0834]
		dischargeSOC0_2C = [0, 0.0864, 0.2006, 0.4012, 0.5000, 0.7377, 0.8519, 0.8951, 0.9259, 0.9475, 0.9630, 0.9784, 0.9907, 1.0000]
		
		chargeVoltage = [3.5450, 3.8018, 3.8404, 3.9046, 3.9495, 3.9881, 4.0459, 4.1358, 4.1936, 4.2000]
		chargeVoltageH = [0, 0.0953, 0.2187, 0.4935, 0.7626, 0.9981, 1.2505, 1.4972, 1.6486, 3.0000]
		chargeCurrent = [1.1000, 1.1000, 0.7167, 0.5833, 0.4333, 0.2667, 0.1500, 0.0667, 0.0333]
		chargeCurrentH = [0, 1.7069, 1.8017, 1.8793, 1.9828, 2.1638, 2.3879, 2.7328, 3.0000]
		chargeCapacity = [0.2333, 0.3500, 1.9000, 1.9500, 2.0667, 2.1500, 2.2333, 2.2667, 2.2667, 2.2833]
		chargeCapacityH = [0, 0.0776, 1.6034, 1.7069, 1.7931, 2.0000, 2.2845, 2.5603, 2.7500, 3.0000]
	end
	methods (Static)
		function plotDischarge()
			figure
			hold on
			plot(LIR18650.dischargeSOC1C, LIR18650.dischargeVoltage1C, '-r');
			plot(LIR18650.dischargeSOC0_5C, LIR18650.dischargeVoltage0_5C, '-g');
			plot(LIR18650.dischargeSOC0_2C, LIR18650.dischargeVoltage0_2C, '-b');
			title ('Typical Discharge Chararcteristics')
			legend ('1C', '0.5C', '0.2C')
			xlabel('SoC [0-1]')
			ylabel('Voltage (V)')
			ylim ([2.6, 4.2])
		end
		
		function plotCharge()
			figure
			hold on
			yyaxis left
			plot(LIR18650.chargeVoltageH, LIR18650.chargeVoltage, '-r');
			ylabel('Voltage (V)')
			ylim ([2.8, 4.3])
			
			yyaxis right
			plot(LIR18650.chargeCurrentH, LIR18650.chargeCurrent, '-b');
			plot(LIR18650.chargeCapacityH, LIR18650.chargeCapacity, '-g');
			ylabel('Current/Charge (A/Ah)')
			ylim ([0, 2.5])

			title ('Typical Charge Chararcteristics')
			legend ('Voltage', 'Current', 'Capacity')
			xlabel('Time (h)')
			
		end
	end
end
