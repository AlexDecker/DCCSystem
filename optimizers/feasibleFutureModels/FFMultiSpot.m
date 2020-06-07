%See FeasibleFuture class for more information. Based on MultiSpot Algorithm

classdef FFMultiSpot < FeasibleFuture
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
        nt
		nr
        cloud
    end
    methods
        function obj = FFMultiSpot(hashSize, nSegments, maxSize, nt, nr)

            obj.hashSize = hashSize;
            obj.nSegments = nSegments;
            obj.maxSize = maxSize;
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
            
            %get any element
            [Q0,D0] = initialSet.cloud.any();

            Rl = FFUtils.calculateLoadResistances(Q0, deviceData, chargeData);
                
            %calculating the minimal receiving current to keep alive
            minIr = FFUtils.calculateMinIr(Q0, deviceData, chargeData,...
                timeSlot.Id, dt);

			if FFMultiSpot.verbose_top
				print_vector('Q0', Q0, 0);
				print_vector('RL', Rl, 0);
				print_vector('minIr', minIr, 0);
			end

            %inverse of the impedance matrix
			Z = timeSlot.Z+diag([zeros(obj.nt,1);Rl]);
			
			%notable sub-matrices
			ZR = Z(obj.nt+1:end,obj.nt+1:end);
			RR = diag(diag(real(ZR)));
			%ir = H*it
			%0 = Z(nt+1:end,:)i
			%0 = Z(nt+1:end,1:nt)it + Z(nt+1:end,nt+1:end)ir
			%Z(nt+1:end,nt+1:end)ir = -Z(nt+1:end,1:nt)it
			%ir = -(Z(nt+1:end,nt+1:end)\Z(nt+1:end,1:nt))*it
			H = -(Z(obj.nt+1:end,obj.nt+1:end)\Z(obj.nt+1:end,1:obj.nt));
			
				
			[Vectors,Values] = eig(H'*RR*H);
			%the receiving power is maximized when the transmitting currents are proportional to the 
			%eigenvector of greatest eigenvalue of H'*RR*H and all power is spent
			[~,i] = max(diag(Values));
			it_base = Vectors(:,i);
			
			ir_base = H*it_base;
			i_base = [it_base;ir_base];
			v_base = Z*i_base;
			v_base = v_base(1:obj.nt);

			%range of voltage multipliers which lead to feasible states
			[minK, maxK] = FFUtils.calculateLimitConstants(v_base,i_base,...
				minIr, constraints);
			
			if FFMultiSpot.verbose
				print_vector('v_base', v_base, 1);
				print_vector('it_base', abs(it_base), 1);
				print_vector('ir_base', abs(ir_base), 1);
				disp(['...', num2str(minK), '<=k<=', num2str(maxK)]);
			end
			
			if minK + FFMultiSpot.tolerance <= maxK - FFMultiSpot.tolerance
				coeff = rand;
				K = coeff * minK + (1-coeff) * maxK;
				%-K is not required to be tested, since it leads to the same abs(ir)

				V = K*v_base;
				I = K*i_base;

				Q = FFUtils.integrateCharge(Q0,I,timeSlot.Id,deviceData,dt,chargeData);
			
				if FFMultiSpot.verbose_down
					print_vector('K', K, 2);
					print_vector('V', V, 2);
					print_vector('It', abs(I(1:obj.nt)), 2);
					print_vector('Ir', abs(I(obj.nt+1:end)), 2);
					print_vector('Q', Q, 2);
					print_vector('Id', timeSlot.Id, 2);
				end
	
				cloud = cloud.insert(Q,V,D0);
			
				%is this element a valid final state?
				if mean(Q >= chargeData.threshold)==1
					final = struct('voltage',V,'previous',D0,'charge',Q);
					if stop_if_threshold_reached
						threshold_reached = true;
					end
				end
            end
            %disp(['New feasible future with ',num2str(cloud.countElements()), ' elements.']);

            %build the object
            new = FFMultiSpot(obj.hashSize, obj.nSegments, obj.maxSize, obj.nt, obj.nr);

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
            initial = FFMultiSpot(1, obj.nSegments, 1, obj.nt, obj.nr);

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
