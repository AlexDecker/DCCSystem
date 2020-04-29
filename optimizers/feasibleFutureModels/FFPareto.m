%See FeasibleFuture class for more information.
%This class implements the feasible future as a large set of elements, each one
%being equivalent to a pareto-optimal reachable state given an initial set of
%states and a given time slot.
classdef FFPareto < FeasibleFuture
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
        thr_top
        thr
        ttl_top
        ttl
        nt
		nr
        cloud
    end
    methods
        function obj = FFPareto(hashSize, nSegments, maxSize, thr_top, thr, ttl_top,...
            ttl, nt, nr)

            obj.hashSize = hashSize;
            obj.nSegments = nSegments;
            obj.maxSize = maxSize;
            obj.thr_top = thr_top;
            obj.thr = thr;
            obj.ttl_top = ttl_top;
            obj.ttl = ttl;
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
            
            %The seach is divided into three hierarquical levels:
            %Top: for each state in initialSet
            %regular (no suffix): for each basic transmitting voltage vector

            %each search stops:
            %-if all recent attempts of generating a new state have failed
            %-if the cloud already has enough content
            %-if the generated sample is too large (it avoid the lack of variability regarding
            %the superior levels. 
            
            %failure: no vector was inserted 
            consecutive_failures_top = 0;
            successes_top = 0;
            attempts_top = 0;
            threshold_reached = false;
            while consecutive_failures_top < obj.thr_top && successes_top < obj.maxSize && ...
                attempts_top < obj.ttl_top && ~threshold_reached && ~cloud.isFull()

                attempts_top = attempts_top + 1;
                
                %get any element
                [Q0,D0] = initialSet.cloud.any();

                Rl = FFUtils.calculateLoadResistances(Q0, deviceData, chargeData);
                
                %calculating the minimal receiving current to keep alive
                minIr = FFUtils.calculateMinIr(Q0, deviceData, chargeData,...
                    timeSlot.Id, dt);

                if FFPareto.verbose_top
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
                while consecutive_failures < obj.thr && successes_top < obj.maxSize &&...
                    attempts < obj.ttl && ~threshold_reached && ~cloud.isFull()
                    
                    attempts = attempts + 1;

                    %the base voltage (some positive, some negative)
                    v_base = rand(obj.nt,1)-0.5;
                    %the base current
                    i_base = iZ*[v_base; zeros(obj.nr,1)];

                    %range of voltage multipliers which lead to feasible states
                    [minK, maxK] = FFUtils.calculateLimitConstants(v_base,i_base,...
                        minIr, constraints);
                    
                    if FFPareto.verbose
						print_vector('v_base', v_base, 1);
						print_vector('it_base', abs(i_base(1:obj.nt)), 1);
						print_vector('ir_base', abs(i_base(obj.nt+1:end)), 1);
                        disp(['...', num2str(minK), '<=k<=', num2str(maxK)]);
                    end

                    if minK + FFPareto.tolerance > maxK - FFPareto.tolerance
                        %the range is empty
                        consecutive_failures = consecutive_failures+1;
                    else
						%-K is not required to be tested, since it leads to the same abs(ir)

						V = (maxK - FFPareto.tolerance)*v_base;
						I = (maxK - FFPareto.tolerance)*i_base;
						
						if max(abs(I) > constraints.maxCurr)==1
							disp('abs(I) should not be over the maximum!');
							continue;
						end
						
						if real(I(1:obj.nt)'*V) > constraints.maxPact
							disp('P should not be over the maximum!');
							continue;
						end

						Q = FFUtils.integrateCharge(Q0,I,timeSlot.Id,deviceData,dt,chargeData);
						
						if FFPareto.verbose_down
							print_vector('K', maxK - FFPareto.tolerance, 2);
							print_vector('V', V, 2);
							print_vector('It', abs(I(1:obj.nt)), 2);
							print_vector('Ir', abs(I(obj.nt+1:end)), 2);
							print_vector('Q', Q, 2);
							print_vector('Id', timeSlot.Id, 2);
						end
						
						if max(Q<=chargeData.minimum)==1
							disp('Q should not be under the minimum!');
							continue;
						end

						%is it already in the cloud?
						D = cloud.discretize(Q);
						[found,~,~,~] = cloud.search(D);
						if found
							%yes, it is
							consecutive_failures = consecutive_failures + 1;
							if FFPareto.verbose_down
								disp('The hyper-cube is already full');
							end
						else
							successes = successes+1;
							consecutive_failures = 0;
							
							cloud = cloud.insert(Q,V,D0);
							
							%is this element a valid final state?
							if mean(Q >= chargeData.threshold)==1
								final = struct('voltage',V,'previous',D0,'charge',Q);
								if stop_if_threshold_reached
									threshold_reached = true;
								end
							end
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
            new = FFPareto(obj.hashSize, obj.nSegments, obj.maxSize, obj.thr_top, obj.thr,...
                obj.ttl_top, obj.ttl, obj.nt, obj.nr);

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
            initial = FFPareto(1, obj.nSegments, 1, obj.thr_top, obj.thr, obj.ttl_top,...
				obj.ttl, obj.nt, obj.nr);

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
end
