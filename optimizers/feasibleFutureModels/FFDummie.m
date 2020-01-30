%See FeasibleFuture class for more information.
%This class implements the feasible future as a large set of elements, each one
%being equivalent to a reachable state given an initial set of states and a given
%time slot.
classdef FFDummie < FeasibleFuture
    properties
        hashSize
        nSegments
        maxSize
        thr_top
        thr
        nt
		nr
        cloud
    end
    methods
        function obj = FFDummie(hashSize,nSegments,maxSize,thr_top,thr,nt,nr)
            obj.hashSize = hashSize;
            obj.nSegments = nSegments;
            obj.maxSize = maxSize;
            obj.thr_top = thr_top;
            obj.thr = thr;
            obj.nt = nt;
			obj.nr = nr;
            obj.cloud = [];
        end
        
        function new = newFeasibleFuture(obj, initialSet, timeSlot, dt,...
            chargeData, deviceData, constraints)
            
            %creating the cloud to store the generated points
            cloud = CloudHash(obj.hashSize, obj.nSegments, chargeData.minimum,...
                chargeData.maximum, obj.maxSize, obj.nt);

            while true
                %get any element
                Q0 = initialSet.any();

                Rl = FFDummie.calculateLoadResistances(Q0, deviceData, chargeData);
                
                %calculating the minimal receiving current to keep alive
                minIr = FFDummie.calculateMinIr(Q0, deviceData, chargeData,...
                    timeSlot.Id, dt);

                %inverse of the impedance matrix
                iZ = eye(obj.nt+obj.nr)\(timeSlot.Z+diag(Rl));
                
                %number of failures: number of times when no vector was inserted due
                %it has already been inserted or minK>maxK
                consecutive_failures_top = 0;
                successes_top = 0;
                %generate a set of new states parting from Q
                while consecutive_failures_top < obj.thr_top && successes_top < obj.maxSize
                    %the base voltage
                    v_base = rand(obj.nt,1);
                    %the base current
                    i_base = iZ*v_base;
                    
                    [minK, maxK] = FFDummie.calculateLimitConstants(v_base,i_base,...
                        minIr, constraints);
                    if minK>maxK
                        consecutive_failures_top = consecutive_failures_top+1;
                    else
                        successes = 0;
                        consecutive_failures = 0;
                        while consecutive_failures < obj.thr && successes_top < obj.maxSize
                            %create a new future state
                            K = rand*(maxK-minK) + minK;
                            V = K*v_base;
                            Q = FFDummie.integrateCharge(Q0,I,Id,deviceData,dt);
                            %is it already in the cloud?
                            D = cloud.discretize(Q);
                            [found,~,~,~] = cloud.search(D);
                            if found
                                %yes, it is
                                consecutive_failures = consecutive_failures+1;
                            else
                                %no, so insert it
                                successes = successes+1;
                                successes_top = successes_top+1;
                                consecutive_failures = 0;
                                cloud = cloud.insert(D,V,D0);
                            end
                        end
                        
                        if successes>0
                            consecutive_failures_top = 0;
                        else
                            consecutive_failures_top = consecutive_failures_top + 1;
                        end
                    end
                end
            end
        end
        
        function initial = generateInitialSet(obj, chargeData)
            %creating the cloud to store the point
            cloud = CloudHash(1, obj.nSegments, chargeData.minimum,...
                chargeData.maximum, 1, obj.nt);

            %insert the initial state
            d = cloud.discretize(chargeData.initial);
            cloud = cloud.insert(d,zeros(obj.nr,1),zeros(obj.nr,1));

            %build the object
            initial = FFDummie(1, obj.nSegments, 1, obj.nt, obj.nr)
        end
        
        %search a given charge vector in the set, returning a structure containing
        %the following fields:
        %   * charge: the chargeVector itself
        %   * voltages: the active voltage vector to turn previous into q
        %   * previous: the charge vector from the initial set
        function q = search(obj, chargeVector)
            q = [];
        end

        %returns true if there is no element in the set
        function b = isEmpty(obj)
            b = true;
        end
    end
    methods(Static)
        function Rl = calculateLoadResistances(Q, deviceData, chargeData)
            %calculating the load resistance of each receiving device
			Rl = zeros(obj.nr);
			for r = 1:obj.nr
				Rl(r) = obj.deviceData(r).getRLfromSOC(q(r)/...
					obj.chargeData.maximum(r));
			end 
        end

        function minIr = calculateMinIr(Q0, deviceData, chargeData, Ir, dt)
            %minimal receiving amplitude for each device to stay alive (exclusive)
            minIr = zeros(obj.nr,1);
            for r=1:obj.nr
                %minimal charge current
                ic = (chargeData.minimum-Q0)/dt;
                %minimal input current (after powering the device)
                input = deviceData(r).iEffectiveChargeCurrent(ic);
                %minimal receiving current amplitude
                minIr = deviceData(r).iConvACDC(input+id);
            end
        end

        function [minK, maxK] = calculateLimitConstants(v_base,i_base,...
            minIr, constraints)
            %i_base: only transmitting currents
            it_base = i_base(1:obj.nt);
            %i_base: only receiving currents
            ir_base = i_base(obj.nt+1:end);
            
            %the maximum k which respects the apparent-power constraint
            maxK = sqrt(constraints.maxPapp/abs((it_base')*v_base));
            %the maximum k which respects the active-power constraint
            maxK = min(maxK,sqrt(constraints.maxPact/real((it_base')*v_base)));

            for i=1:obj.nt+obj.nr
                %the maximum k which respects the current constraint for element i
                maxK = min(maxK,constraints.maxCurr(i)/abs(i_base(i)));
            end
            
            minK = inf;
            for r=1:obj.nr
                %the minimum k which respects the minIr (exclusive)
                minK = max(minK, minIr/abs(ir_base(r)));
            end
        end

        function Q = integrateCharge(Q0, I, Id, deviceData, dt)
            %receiving amplitudes
            Ir = abs(I(obj.nt+1:end));    
            
            Q = zeros(obj.nr,1);
            for r=1:obj.nr
                %received current (CC)
                ir = deviceData(r).convACDC(Ir(r));
                %charging input current
                ic = ir-Id(r);
                %final charge
                Q(r) = deviceData(r).effectiveChargeCurrent(ic)*dt+Q0(r); 
            end
        end
    end
end
