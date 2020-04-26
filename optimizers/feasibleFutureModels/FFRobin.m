%See FeasibleFuture class for more information.
%This class implements the feasible future as a unitary set of elements. The element is
%the one which maximizes the minimum charge. If stop_if_threshold_reached=true, it
%also searches for a final solution
classdef FFRobin < FeasibleFuture
    properties(Constant)
        tolerance = 1e-6
        verbose_top = false
        verbose = false
        verbose_down = false
    end
    properties
        hashSize
        nSegments
        maxSize
        thr
        thr_down
        ttl
        ttl_down
        nt
		nr
        cloud
    end
    methods
        function obj = FFRobin(hashSize, nSegments, maxSize, thr, thr_down,...
            ttl, ttl_down, nt, nr)

            obj.hashSize = hashSize;
            obj.nSegments = nSegments;
            obj.maxSize = maxSize;
            obj.thr = thr;
            obj.thr_down = thr_down;
            obj.ttl = ttl;
            obj.ttl_down = ttl_down;
            obj.nt = nt;
			obj.nr = nr;
            obj.cloud = [];
			
        end
        
        function [final, new] = newFeasibleFuture(obj, initialSet, timeSlot, dt,...
            chargeData, deviceData, constraints, stop_if_threshold_reached)
            
            final = [];

            %creating the cloud to store the generated points
            cloud = DoubleCloudHash(obj.hashSize, obj.nSegments, chargeData.minimum,...
                chargeData.maximum, obj.maxSize, obj.nt);
				
			best_Q = -ones(obj.nr,1);%default: invalid
			best_D0 = zeros(obj.nr,1);
			best_V = zeros(obj.nt,1);
            
            %The seach is divided into three hierarquical levels:
            %regular (no suffix): for each basic transmitting voltage vector
            %Down: for each multiple of the basic transmitting voltage vector

            %each search stops:
            %-if all recent attempts of generating a new state have failed
            %-if the cloud already has enough content
            %-if the generated sample is too large (it avoid the lack of variability regarding
            %the superior levels. 
            
            %failure: no vector was inserted 
                
			%get any element
			[Q0,D0] = initialSet.cloud.any();

			Rl = FFRobin.calculateLoadResistances(Q0, deviceData, chargeData);
			
			%calculating the minimal receiving current to keep alive
			minIr = FFRobin.calculateMinIr(Q0, deviceData, chargeData,...
				timeSlot.Id, dt);

			if FFRobin.verbose_top
				print_vector('Q0', Q0, 0);
				print_vector('RL', Rl, 0);
				print_vector('minIr', minIr, 0);
			end

			%inverse of the impedance matrix
			Z = timeSlot.Z+diag([zeros(obj.nt,1);Rl]);
			iZ = eye(obj.nt+obj.nr)/Z;
			
			%number of failures: number of times when no vector was inserted due
			%it has already been inserted or minK>maxK
			consecutive_failures = 0;
			successes = 0;
			attempts = 0;
			%generate a set of new states parting from Q
			while consecutive_failures < obj.thr && attempts < obj.ttl &&...
				~threshold_reached && ~cloud.isFull()
				
				attempts = attempts + 1;

				%the base voltage (some positive, some negative)
				v_base = rand(obj.nt,1)-0.5;
				%the base current
				i_base = iZ*[v_base; zeros(obj.nr,1)];

				%range of voltage multipliers which lead to feasible states
				[minK, maxK] = FFRobin.calculateLimitConstants(v_base,i_base,...
					minIr, constraints);
				
				if FFRobin.verbose
					print_vector('v_base', v_base, 1);
					print_vector('it_base', abs(i_base(1:obj.nt)), 1);
					print_vector('ir_base', abs(i_base(obj.nt+1:end)), 1);
					disp(['...', num2str(minK), '<=k<=', num2str(maxK)]);
				end

				if minK + FFRobin.tolerance > maxK - FFRobin.tolerance
					%the range is empty
					consecutive_failures = consecutive_failures+1;
				else
					successes_down = 0;
					consecutive_failures_down = 0;

					for K = linspace(max(maxK - FFRobin.tolerance, minK + FFRobin.tolerance),...
						min(maxK - FFRobin.tolerance, minK + FFRobin.tolerance), obj.ttl_down)
						
						if consecutive_failures_down >= obj.thr_down || threshold_reached || cloud.isFull()
							break;
						end

						%-K is not required to be tested, since it leads to the same abs(ir)

						V = K*v_base;
						I = K*i_base;

						Q = FFRobin.integrateCharge(Q0,I,timeSlot.Id,deviceData,dt,chargeData);
						
						if FFRobin.verbose_down
							print_vector('K', K, 2);
							print_vector('V', V, 2);
							print_vector('It', abs(I(1:obj.nt)), 2);
							print_vector('Ir', abs(I(obj.nt+1:end)), 2);
							print_vector('Q', Q, 2);
							print_vector('Id', timeSlot.Id, 2);
						end

						successes_down = successes_down+1;
						successes = successes+1;
						consecutive_failures_down = 0;
						
						if min(Q)>min(best_Q)
							best_D0 = D0;
							best_V = V;
						end
						
						%is this element a valid final state?
						if stop_if_threshold_reached
							if mean(Q >= chargeData.threshold)==1
								final = struct('voltage',V,'previous',D0,'charge',Q);
								threshold_reached = true;
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
			
			if threshold_reached
				%there is a final solution
				cloud = cloud.insert(final.Q,final.V,final.D0);
			else
				%if there is at least a valid solution
				if min(Q_best>=0)==1
					%insert it
					cloud = cloud.insert(Q_best,V_best,D0_best);
				end
			end
            
            disp(['New feasible future with ',num2str(cloud.countElements()), ' elements.']);

            %build the object
            new = FFRobin(obj.hashSize, obj.nSegments, obj.maxSize, obj.thr,...
                obj.thr_down, obj.ttl, obj.ttl_down, obj.nt, obj.nr);

            %insert the cloud into the object
            new.cloud = cloud;
        end
        
        function initial = generateInitialSet(obj, chargeData)
            %creating the cloud to store the point
            cloud = DoubleCloudHash(1, obj.nSegments, chargeData.minimum,...
                chargeData.maximum, 1, obj.nt);

            %insert the initial state
            cloud = cloud.insert(chargeData.initial, zeros(obj.nt,1), zeros(obj.nr,1));

            %build the object
            initial = FFRobin(1, obj.nSegments, 1, obj.thr, obj.thr_down,...
                obj.ttl, obj.ttl_down, obj.nt, obj.nr);

            %insert the cloud into the object
            initial.cloud = cloud;
        end
        
        %search a given discretized charge vector in the set, returning a
        %structure containing the following fields:
        %   * voltage: the active voltage vector to turn previous into q
        %   * previous: the discretized charge vector from the initial set
        function element = search(obj, dChargeVector)
            [found, h, j, pj] = obj.cloud.search(dChargeVector);
            if found
                if isempty(j)
                    [~,Q,D0,V] = obj.cloud.readFromPool(pj);
                else
                    [~,Q,D0,V] = obj.cloud.read(h,j);
                end
                element = struct('voltage',V,'previous',D0,'charge',Q);
            else
                element = [];
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
				if SOC < 0
					error('SOC must be no less than 0');
				elseif SOC > 1
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
                            if input < minIn
                                %trivial, since the output of the conversion is non-negative
                                minIr(r) = -inf;
                            else
                                %minimal receiving current amplitude
                                minIr(r) = deviceData(r).iConvACDC(input);
                            end
                        end
                    end
                end
            end
            %transforming the limit from exclusive to inclusive by adding a very small
            %tolerance value
            minIr = minIr + FFRobin.tolerance;
            %avoid negative "amplitudes"
            minIr = max(0, minIr);
        end

        function [minK, maxK] = calculateLimitConstants(v_base,i_base,...
            minIr, constraints)
            %i_base: only transmitting currents
            it_base = i_base(1:length(v_base));
            %i_base: only receiving currents
            ir_base = i_base(length(v_base)+1:end);
            
            %the maximum k which respects the active-power constraint
            maxK = sqrt((constraints.maxPact - FFRobin.tolerance)/real((it_base')*v_base));

            for i=1:length(i_base)
                %the maximum k which respects the current constraint for element i
                maxK = min(maxK,(constraints.maxCurr(i) - FFRobin.tolerance)/abs(i_base(i)));
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
                Q(r) = min(ic(r)*dt+Q0(r), chargeData.maximum(r)-FFRobin.tolerance);
				Q(r) = max(0, Q(r));
            end
			
			if FFRobin.verbose_down
				disp('Conversions: begin');
				print_vector('Receiving current', Ir, 2);
				print_vector('Input current', ir, 2);
				print_vector('Charging current', ic, 2);
				disp('Conversions: end');
			end
        end
    end
end
