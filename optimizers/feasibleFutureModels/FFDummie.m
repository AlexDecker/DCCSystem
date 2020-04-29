%See FeasibleFuture class for more information.
%This class implements the feasible future as a large set of elements, each one
%being equivalent to a reachable state given an initial set of states and a given
%time slot.
classdef FFDummie < FeasibleFuture
    properties(Constant)
        tolerance = 1e-6
		tolerance_fine_adjustment = 1e-9
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
        thr_down
        ttl_top
        ttl
        ttl_down
		ttl_adjustment
        nt
		nr
        cloud
    end
    methods
        function obj = FFDummie(hashSize, nSegments, maxSize, thr_top, thr, thr_down, ttl_top,...
            ttl, ttl_down, ttl_adjustment, nt, nr)

            obj.hashSize = hashSize;
            obj.nSegments = nSegments;
            obj.maxSize = maxSize;
            obj.thr_top = thr_top;
            obj.thr = thr;
            obj.thr_down = thr_down;
            obj.ttl_top = ttl_top;
            obj.ttl = ttl;
            obj.ttl_down = ttl_down;
			obj.ttl_adjustment = ttl_adjustment;
            obj.nt = nt;
			obj.nr = nr;
            obj.cloud = [];
			
        end
        
        function [final, new] = newFeasibleFuture(obj, initialSet, timeSlot, dt,...
            chargeData, deviceData, constraints, stop_if_threshold_reached)
            
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

                if FFDummie.verbose_top
					disp(['Attempt number ',num2str(attempts_top)]);
					disp(['Successes (top) ',num2str(successes_top)]);
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
                    
                    if FFDummie.verbose
						print_vector('v_base', v_base, 1);
						print_vector('it_base', abs(i_base(1:obj.nt)), 1);
						print_vector('ir_base', abs(i_base(obj.nt+1:end)), 1);
                        disp(['...', num2str(minK), '<=k<=', num2str(maxK)]);
                    end

                    if minK + FFDummie.tolerance > maxK - FFDummie.tolerance
                        %the range is empty
                        consecutive_failures = consecutive_failures+1;
                    else
                        successes_down = 0;
                        consecutive_failures_down = 0;

						for K = linspace(max(maxK - FFDummie.tolerance, minK + FFDummie.tolerance),...
							min(maxK - FFDummie.tolerance, minK + FFDummie.tolerance), obj.ttl_down)
							
							if consecutive_failures_down >= obj.thr_down ||...
								successes_top >= obj.maxSize || threshold_reached || cloud.isFull()
								break;
							end

                            %-K is not required to be tested, since it leads to the same abs(ir)

                            V = K*v_base;
                            I = K*i_base;

                            Q = FFUtils.integrateCharge(Q0,I,timeSlot.Id,deviceData,dt,chargeData);
                            
                            if FFDummie.verbose_down
								print_vector('K', K, 2);
								print_vector('V', V, 2);
								print_vector('It', abs(I(1:obj.nt)), 2);
								print_vector('Ir', abs(I(obj.nt+1:end)), 2);
								print_vector('Q', Q, 2);
								print_vector('Id', timeSlot.Id, 2);
                            end

                            %is it already in the cloud?
                            D = cloud.discretize(Q);
                            [found,~,~,~] = cloud.search(D);
                            if found
                                %yes, it is
                                consecutive_failures_down = consecutive_failures_down+1;
								if FFDummie.verbose_down
									disp('The hyper-cube is already full');
								end
                            else
                                %The next slot will consider only the center of the hypercube defined
								%by D. Q belongs to hypercube D, but the difference between the two 
								%stated may lead to errors. In particular, the error increases linearly
								%regarding the time. So, we adjust the voltage so the next state becames
								%the center of D
								
								[minQ, maxQ] = cloud.dediscretize(D);
								Q_center = (minQ + maxQ)/2; %the center of D
								
								if FFDummie.verbose_down
									print_vector('Q_center',Q_center,2);
								end
								
								%In first place: is this state feasible regarding minimum charge?
								if min(Q_center > chargeData.minimum)==0
									%no, so it is a failure
									consecutive_failures_down = consecutive_failures_down+1;
									if FFDummie.verbose_down
										disp('The center of the hyper-cube is not feasible due the minimum charge constraint.');
									end
								else
									Ic_center = (Q_center - Q0)/dt; %the required effective charge current
									if FFDummie.verbose_down
										print_vector('Ic_center', Ic_center, 2);
									end
									target_reacheable = true; %default
									%some functions have constant intervals in their domain, so their inverse is not a function.
									%However, as they are monotonically increasing, we can manage the inverse as a function whose
									%image is an interval.
									min_targetIr = zeros(obj.nr,1);
									max_targetIr = zeros(obj.nr,1);
									for r = 1:obj.nr
										[min_Ic, max_Ic] = deviceData(r).domain_iEffectiveChargeCurrent();
										%is this effective charge current possible?
										if min_Ic <= Ic_center(r) && Ic_center(r) <= max_Ic
											%the required input DC current
											[In0,In1] = deviceData(r).iEffectiveChargeCurrent(Ic_center(r));
											In0 = In0 + timeSlot.Id(r) - FFDummie.tolerance;
											In1 = In1 + timeSlot.Id(r) + FFDummie.tolerance;
											
											%the interest region of the domain
											[min_In, max_In] = deviceData(r).domain_iConvACDC();
											min_In = max(min_In, In0);
											max_In = min(max_In, In1);
											if min_In <= max_In
												%the required amplitude for the receiving current
												[min_targetIr(r), ~] = deviceData(r).iConvACDC(min_In);
												
												%verifying if the maximum current constraint is satisfied
												if min_targetIr(r) + FFDummie.tolerance > constraints.maxCurr(obj.nt+r)
													target_reacheable = false;
												else
													[~, max_ir] = deviceData(r).iConvACDC(max_In);
													max_targetIr(r) = min(max_ir, constraints.maxCurr(obj.nt+r) - FFDummie.tolerance);
												end
												
											else
												target_reacheable = false;
												break;
											end
										else
											target_reacheable = false;
											break;
										end
									end
									
									if target_reacheable									
										%now for targetIr itself we must choose the vector inside the interval [min_targetIr, max_target_Ir]
										%which is the closest to the former receiving voltage vector, that is, abs(I(obj.nt+1:end))
										
										%first: which currents already are inside the target interval?
										inside = abs(I(obj.nt+1:end)) <= max_targetIr & abs(I(obj.nt+1:end)) >= min_targetIr;
										%which ones are under the interval?
										under = abs(I(obj.nt+1:end)) < min_targetIr;
										%what about over?
										over = abs(I(obj.nt+1:end)) > max_targetIr;
										
										targetIr = inside.*abs(I(obj.nt+1:end)) + under.*min_targetIr + over.*max_targetIr;
										
										if FFDummie.verbose_down
											print_vector('targetIr', targetIr, 2);
										end
									
										%search for the little adjustment of the voltage which will minimize the errors
										[dv,~] = fine_adjustment(Z, V, I(1:obj.nt), I(obj.nt+1:end),...
											constraints.maxCurr(1:obj.nt) - FFDummie.tolerance, targetIr,...
											constraints.maxPact - FFDummie.tolerance,...
											FFDummie.tolerance_fine_adjustment, obj.ttl_adjustment);
										
										if isempty(dv)
										
											% the returned values are not valid
											consecutive_failures_down = consecutive_failures_down+1;
											
											if FFDummie.verbose_down
												disp('Fine-Adjustment failure!!!');
											end
											
										else
											
											V_new = V + dv;
											
											%just verifying the returned values..
											I_test = iZ*[V_new;zeros(obj.nr,1)];
											It_test = I_test(1:obj.nt);
											P_test = V_new.'*real(It_test);
											
											if sum(abs(It_test)>constraints.maxCurr(1:obj.nt)) > 0 || ...
												P_test > constraints.maxPact || ...
												sum(abs(targetIr - abs(I_test(obj.nt+1:end))) > FFDummie.tolerance) > 0
												
												% the returned values are not valid
												consecutive_failures_down = consecutive_failures_down+1;
												
												if FFDummie.verbose_down
													disp('Fine-Adjustment error');
												end
												
											else
											
												successes_down = successes_down+1;
												successes = successes+1;
												consecutive_failures_down = 0;
												
												D_new = cloud.discretize(Q_center);
												cloud = cloud.insert(D_new,V_new,D0);
												
												%is this element a valid final state?
												if mean(Q_center>=chargeData.threshold)==1
													final = struct('voltage',V_new,'previous',D0,'charge',Q_center);
													if stop_if_threshold_reached
														threshold_reached = true;
													end
												end
											end
										end
									else
										%target is unreacheable. failure
										consecutive_failures_down = consecutive_failures_down+1;
										
										if FFDummie.verbose_down
											disp('The center of the hyper-cube is not feasible.');
										end
									end
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
					successes_top = successes_top + 1;
                else
                    %no successes in the child-level
                    consecutive_failures_top = consecutive_failures_top + 1;
                end
            end
            
            disp(['New feasible future with ',num2str(cloud.countElements()), ' elements.']);

            %build the object
            new = FFDummie(obj.hashSize, obj.nSegments, obj.maxSize, obj.thr_top, obj.thr,...
                obj.thr_down, obj.ttl_top, obj.ttl, obj.ttl_down, obj.ttl_adjustment,...
				obj.nt, obj.nr);

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
                obj.ttl_top, obj.ttl, obj.ttl_down, obj.ttl_adjustment, obj.nt, obj.nr);

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
                    [D,D0,V] = obj.cloud.readFromPool(pj);
                else
                    [D,D0,V] = obj.cloud.read(h,j);
                end
				[minQ,maxQ] = obj.cloud.dediscretize(D);
                d = struct('voltage',V,'previous',D0,'charge',(minQ + maxQ) / 2);
            else
                d = [];
            end
        end

        %returns true if there is no element in the set
        function b = isEmpty(obj)
            b = (obj.cloud.countElements()==0);
        end
    end
end
