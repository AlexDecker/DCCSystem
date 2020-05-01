%See FeasibleFuture class for more information.
%This class implements the feasible future as a large set of elements
classdef FFUtils
    properties(Constant)
        tolerance = 1e-6
    end
    methods
        function obj = FFUtils()			
        end
    end
    methods(Static)
        function Rl = calculateLoadResistances(Q, deviceData, chargeData)
            %calculating the load resistance of each receiving device
			Rl = zeros(length(Q),1);
			for r = 1:length(Q)
                SOC = Q(r)/chargeData.maximum(r);
				if SOC < 0
					SOC
					error('SOC must be no less than 0');
				elseif SOC > 1
					SOC
					error('SOC must be no more than 1');
				end
				Rl(r) = deviceData(r).getRLfromSOC(SOC);
			end 
        end

        function minIr = calculateMinIr(Q0, deviceData, chargeData, Id, dt)
            %minimal receiving amplitude for each device to stay alive (exclusive)
            minIr = zeros(length(Q0),1);
            for r=1:length(Q0)
                
                %maximum input current the system is able to provide
                [minIn, maxIn] = deviceData(r).domain_iConvACDC();

                %maximum charging current the system is able to provide
                %(considering Id as the discharge current)
                maxIc = deviceData(r).effectiveChargeCurrent(maxIn-Id(r));

                %minimal charge current (exclusive)
                ic = (chargeData.minimum(r)-Q0(r))/dt;
                
				type = 0;
                if ic >= maxIc
                    %impossible to provide
                    minIr(r) = inf;
					type = 1;
                else
                    if ic < deviceData(r).effectiveChargeCurrent(-Id(r))
                        %trivial (even if there is no receiving current the device is ok)
                        minIr(r) = -inf;
						type = 2;
                    else
                        %minimal input current
                        [input,~] = deviceData(r).iEffectiveChargeCurrent(ic);
						input = input + Id(r);
                        if input >= maxIn
                            %impossible to provide
                            minIr(r) = inf;
							type = 3;
                        else
                            if input < minIn
                                %trivial, since the output of the conversion is non-negative
                                minIr(r) = -inf;
								type = 4;
                            else
                                %minimal receiving current amplitude
                                [minIr(r),~] = deviceData(r).iConvACDC(input);
								type = 5;
                            end
                        end
                    end
                end
				%transforming the limit from exclusive to inclusive by adding a very small
				%tolerance value
				minIr(r) = minIr(r) + FFUtils.tolerance;
				%avoid negative "amplitudes"
				minIr(r) = max(0, minIr(r));
            end
        end

        function [minK, maxK] = calculateLimitConstants(v_base,i_base,...
            minIr, constraints)
            %i_base: only transmitting currents
            it_base = i_base(1:length(v_base));
            %i_base: only receiving currents
            ir_base = i_base(length(v_base)+1:end);
            
            %the maximum k which respects the active-power constraint
            maxK = sqrt((constraints.maxPact - FFUtils.tolerance)/real((it_base')*v_base));

            for i=1:length(i_base)
                %the maximum k which respects the current constraint for element i
                maxK = min(maxK,(constraints.maxCurr(i) - FFUtils.tolerance)/abs(i_base(i)));
            end
            
            minK = -inf;
            for r=1:length(ir_base)
                %the minimum k which respects the minIr (exclusive)
                minK = max(minK, minIr(r)/abs(ir_base(r)));
            end
        end

        function Q = integrateCharge(Q0, I, Id, deviceData, dt, chargeData)
            %receiving amplitudes
            Ir = abs(I(end-length(Id)+1:end));    
            
            Q = zeros(length(Id),1);
			ir = zeros(length(Id),1);
			ic = zeros(length(Id),1);
			
            for r=1:length(Id)
                %received current (CC)
                ir(r) = deviceData(r).convACDC(Ir(r));
				
                %charging input current
                ic(r) = deviceData(r).effectiveChargeCurrent(ir(r)-Id(r));
				
                %final charge
                Q(r) = min(ic(r)*dt+Q0(r), chargeData.maximum(r));%-FFUtils.tolerance);
				Q(r) = max(0, Q(r));
            end
        end
    end
end
