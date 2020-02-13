%See FeasibleFuture class for more information.
%This class implements the feasible future as a large set of elements, each one
%being equivalent to a reachable state given an initial set of states and a given
%time slot.
classdef FFDummie < FeasibleFuture
    properties(Constant)
        tolerance = 1e-6
    end
    properties
        hashSize
        nSegments
        maxSize
        thr_top
        thr
        thr_down
        ttl
        ttl_down
        nt
		nr
        cloud
    end
    methods
        function obj = FFDummie(hashSize,nSegments,maxSize,thr_top,thr,thr_down,ttl,ttl_down,nt,nr)
            obj.hashSize = hashSize;
            obj.nSegments = nSegments;
            obj.maxSize = maxSize;
            obj.thr_top = thr_top;
            obj.thr = thr;
            obj.thr_down = thr_down;
            obj.ttl = ttl;
            obj.ttl_down = ttl_down;
            obj.nt = nt;
			obj.nr = nr;
            obj.cloud = [];
        end
        
        function [final, new] = newFeasibleFuture(obj, initialSet, timeSlot, dt,...
            chargeData, deviceData, constraints)
            
            final = [];

            %creating the cloud to store the generated points
            cloud = CloudHash(obj.hashSize, obj.nSegments, chargeData.minimum,...
                chargeData.maximum, obj.maxSize, obj.nt);
            
            %The seach is divided into three hierarquical levels:
            %Top: for each state in initialSet
            %regular (no suffix): for each basic transmitting voltage vector
            %Down: for each multiple of the basic transmitting voltage vector

            %each search stops:
            %-if all recent attempts of generating a new state have failed
            %-if the cloud already has enough content
            %-if the generated sample is too large (it avoid the lack of variability regarding
            %the superior levels. 
            
            %failure: no vector was inserted 
            consecutive_failures_top = 0;
            successes_top = 0;
            while consecutive_failures_top < obj.thr_top && successes_top < obj.maxSize
                %get any element
                Q0 = initialSet.cloud.any();

                Rl = FFDummie.calculateLoadResistances(Q0, deviceData, chargeData);
                
                %calculating the minimal receiving current to keep alive
                minIr = FFDummie.calculateMinIr(Q0, deviceData, chargeData,...
                    timeSlot.Id, dt);

                %inverse of the impedance matrix
                iZ = eye(obj.nt+obj.nr)\(timeSlot.Z+diag([zeros(obj.nt,1);Rl]));
                
                %number of failures: number of times when no vector was inserted due
                %it has already been inserted or minK>maxK
                consecutive_failures = 0;
                successes = 0;
                attempts = 0;
                %generate a set of new states parting from Q
                while consecutive_failures <i obj.thr && successes_top < obj.maxSize && attempts < obj.ttl
                    
                    attempts = attempts + 1;

                    %the base voltage
                    v_base = rand(obj.nt,1);
                    %the base current
                    i_base = iZ*v_base;
                    
                    %range of voltage multipliers which lead to feasible states
                    [minK, maxK] = FFDummie.calculateLimitConstants(v_base,i_base,...
                        minIr, constraints);

                    if minK>maxK
                        %the range is empty
                        consecutive_failures = consecutive_failures+1;
                    else
                        successes_down = 0;
                        consecutive_failures_down = 0;
                        attempts_down = 0;
                        while consecutive_failures_down < obj.thr_down &&...
                            successes_top < obj.maxSize && attempts_down < obj.ttl_down

                            %create a new future state
                            K = rand*(maxK-minK) + minK;
                            %K*v_base is valid, so is -K*v_base
                            K = sign(rand-0.5)*K;
                            V = K*v_base;
                            Q = FFDummie.integrateCharge(Q0,I,Id,deviceData,dt);
                            %is it already in the cloud?
                            D = cloud.discretize(Q);
                            [found,~,~,~] = cloud.search(D);
                            if found
                                %yes, it is
                                consecutive_failures_down = consecutive_failures_down+1;
                            else
                                %no, so insert it
                                successes_down = successes_down+1;
                                successes = successes+1;
                                consecutive_failures_down = 0;
                                cloud = cloud.insert(D,V,D0);
                                %is this element a valid final state?
                                if mean(Q>=chargeData.threshold)==1
                                    final = struct('voltage',V,'previous',D0);
                                end
                            end
                        end
                        
                        if successes_down > 0
                            %there were some successes in the child-level
                            consecutive_failures = 0;
                        else
                            %no successes in the child-level
                            consecutive_failures = consecutive_failures + 1;
                        end
                    end
                end

                if successes > 0
                    %there were some successes in the child-level
                    consecutive_failures_top = 0;
                else
                    %no successes in the child-level
                    consecutive_failures_top = consecutive_failures_top + 1;
                end
            end
            
            disp(['New feasible future with ',num2str(cloud.countElements()), ' elements.']);

            %build the object
            new = FFDummie(obj.hashSize, obj.nSegments, obj.maxSize, obj.thr_top, obj.thr,...
                obj.thr_down, obj.ttl, obj.ttl_down, obj.nt, obj.nr);

            %insert the cloud into the object
            new.cloud = cloud;
        end
        
        function initial = generateInitialSet(obj, chargeData)
            %creating the cloud to store the point
            cloud = CloudHash(1, obj.nSegments, chargeData.minimum,...
                chargeData.maximum, 1, obj.nt);

            %insert the initial state
            d = cloud.discretize(chargeData.initial);
            cloud = cloud.insert(d,zeros(obj.nt,1),zeros(obj.nr,1));

            %build the object
            initial = FFDummie(1, obj.nSegments, 1, obj.thr_top, obj.thr, obj.thr_down,...
                obj.ttl, obj.ttl_down, obj.nt, obj.nr);

            %insert the cloud into the object
            initial.cloud = cloud;
        end
        
        %search a given discretized charge vector in the set, returning a
        %structure containing the following fields:
        %   * voltage: the active voltage vector to turn previous into q
        %   * previous: the discretized charge vector from the initial set
        function d = search(obj, dChargeVector)
            [found, h, j, pj] = obj.cloud.search(dChargeVector);
            if found
                if isempty(j)
                    [~,D0,V] = obj.cloud.readFromPool(pj);
                else
                    [~,D0,V] = obj.cloud.read(h,j);
                end
                d = struct('voltage',V,'previous',D0)
            else
                d = [];
            end
        end

        %returns true if there is no element in the set
        function b = isEmpty(obj)
            b = (obj.cloud.countElements()==0);
        end
    end
    methods(Static)
        function Rl = calculateLoadResistances(Q, deviceData, chargeData)
            %calculating the load resistance of each receiving device
			Rl = zeros(length(Q),1);
			for r = 1:length(Q)
                SOC = Q(r)/chargeData.maximum(r);
				Rl(r) = deviceData(r).getRLfromSOC(SOC);
			end 
        end

        function minIr = calculateMinIr(Q0, deviceData, chargeData, Id, dt)
            %minimal receiving amplitude for each device to stay alive (exclusive)
            minIr = zeros(length(Q0),1);
            for r=1:length(Q0)
                
                %maximum input current the system is able to provide
                maxIn = deviceData(r).convACDC(deviceData(r).maxReceivingCurrent());

                %maximum charging current the system is able to provide
                %(considering Id as the discharge current)
                maxIc = deviceData(r).effectiveChargeCurrent(maxIn-Id(r));

                %minimal charge current (exclusive)
                ic = (chargeData.minimum(r)-Q0(r))/dt;

                if ic >= maxIc
                    %impossible to provide
                    minIr(r) = inf;
                else
                    if ic < -Id(r)
                        %trivial (even if there is no receiving current the device is ok)
                        minIr(r) = -inf;
                    else
                        %minimal input current
                        input = deviceData(r).iEffectiveChargeCurrent(ic) + Id(r);
                        if input >= maxIn
                            %impossible to provide
                            minIr(r) = inf;
                        else
                            if input < 0
                                %trivial, since the output of the conversion is non-negative
                                minIr(r) = -inf;
                            else
                                %minimal receiving current amplitude
                                minIr(r) = deviceData(r).iConvACDC(input+id);
                            end
                        end
                    end
                end
            end
            %transforming the limit from exclusive to inclusive by adding a very small
            %tolerance value
            minIr = minIr + FFDummie.tolerance;
            %avoid negative "amplitudes"
            minIr = max(0, minIr);
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
