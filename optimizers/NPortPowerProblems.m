%the versions of the problems where the solver is aware about all parameters in the
%system. Minimize the charge time for a given WPT system or define a voltage progression
%for ensuring all devices ramain active for at least maxTau units of time
classdef NPortPowerProblems
    properties
        feasibleFutureModel %used to allow multiple algorithms to manage the 	
		%feasibleFutures

        %time variant data
        timeLine
        dt
        
        %charge constraints and specific data
        chargeData

        %device behavior specific data
        deviceData

        %general constraints (not regarding charge)
        constraints

        %if true, print debug information
        verbose
    end
    properties(SetAccess=private, GetAccess=public)
        n %number of devices
        nt %number of transmitting devices
        nr %number of receiving devices
        maxTau %number of time slots
    end
    properties(Constant)
        plot_colors = 'rgbkmcy'
		tolerance = 1e-7%1e-10
		security_gap = 1e-6
    end
    methods
        %timeLine: data for each time slot. Vector of structs with the following 
		%	fields:
        %   *Z: Impedance matrix without considering saturation or 
		%	charge-dependent resistance
        %   *Id: discharge current (DC, real scalar)
        %dt: duration of one timeslot (real scalar)
        %chargeData: struct with vectors containing data about charge levels
        %   *minimum: vector with the minimum value of charge each device can 
		%	handdle.
        %   *initial: vector with the initial charge of each device
        %   *threshold: vector with the minimum charge of each device at the end
        %   *maximum: the maximum charge supported by each battery
        %deviceData: DeviceData object list (one for each receiver)
        %   *rlCellArray: Lookup table for load resistance given the SOC
        %   *convCellArray: Lookup table for DC current given the AC input 
		%	amplitude, from 0A to maxCurr
        %   *chargeCellArray: Lookup table for the effective charge current givent 
		%	the dc input (from 0A to tha max converted current)
        %constraints: struct contaning
        %   *maxCurr: vector which contains the maximum current amplitude for each 
		%	device
        %   *maxPact: maximum allowed active power
        %feasibleFutureModel: an empty object from a class which inherits from 
		%	FeasibleFuture which is used to create new objects and manage the 
		%	different algorithms according to user's will
        function obj = NPortPowerProblems(timeLine, dt, chargeData, deviceData,...
            constraints, feasibleFutureModel)
            obj.timeLine = timeLine;
            obj.dt = dt;
            obj.chargeData = chargeData;
            obj.deviceData = deviceData;
            obj.constraints = constraints;
            obj.feasibleFutureModel = feasibleFutureModel;
            obj.verbose = false;
            obj = check(obj);
        end
        %check if this object is acceptable
        function obj = check(obj)
            [n1,n2] = size(obj.timeLine);
            if n1==0 || n2~=1
                error('Unexpected dimensions of the timeLine');
            end
            obj.maxTau = n1;%number of time slots
            %verifying the number of devices
            [n1,n2] = size(obj.timeLine(1).Z);
            if n1~=n2
                error('Z matrices must be square');
            end
            obj.n = n1; %total
            [n1,n2] = size(obj.timeLine(1).Id);
            if n1>=obj.n || n2~=1
                error(['Id must be a column vector with one value for each ',...
					'receiving device']);
            end
            obj.nr = n1;%receiving
            obj.nt = obj.n-obj.nr; %transmitting

            %verifying the timeline
            [n1,n2] = size(obj.timeLine);
            if n1==0 || n2~=1
                error('Invalid format for timeLine');
            end
            
            %max discharge current for each device along the timeline
            maxId = zeros(obj.nr,1);
            for t = 1:length(obj.timeLine)
                [n1,n2] = size(obj.timeLine(t).Z);
                if obj.n~=n1 || obj.n~=n2
                    error(['Incompatible size for Z matrix when t=',num2str(t)]);
                end
                [n1,n2] = size(obj.timeLine(t).Id);
                if obj.nr~=n1 || n2~=1
                    error(['Incompatible size for Id when t=',num2str(t)]);
                end
                if sum(sum(real(obj.timeLine(t).Z)<0))>0
                    error(['No negative resistance is allowed. t=',num2str(t)]);
                end
                if sum(sum(real(obj.timeLine(t).Z-diag(diag(...
					obj.timeLine(t).Z)))~=0))>0
                    error(['No real values outside the main diagonal. t=',...
						num2str(t)]);
                end
                %if sum(imag(diag(obj.timeLine(t).Z))~=0)>0
                %    error(['Resonance is required. t=',num2str(t)]);
                %end
                if sum(sum(obj.timeLine(t).Z-obj.timeLine(t).Z.'~=0))>0
                    error(['The mutual inductance must be symmetrical. t=',...
						num2str(t)]);
                end
                if sum(abs(obj.timeLine(t).Id)~=obj.timeLine(t).Id)>0
                    error(['The discharge currents must be real positive. t=',...
						num2str(t)]);
                end
                maxId = max(maxId, obj.timeLine(t).Id);
            end
            
            [n1,n2] = size(obj.dt);
            if n1~=1 || n2~=1 || real(obj.dt)~=obj.dt || obj.dt<=0
                error('dt must be a real positive escalar.');
            end
            
            [n1,n2] = size(obj.chargeData.minimum);
            if n1~=obj.nr || n2~=1 || sum(real(obj.chargeData.minimum)~=...
				obj.chargeData.minimum)>0 || sum(obj.chargeData.minimum<0)>0
                error(['chargeData.minimum must have a real non-negative ',...
					'value for each receiver.']);
            end

            [n1,n2] = size(obj.chargeData.initial);
            if n1~=obj.nr || n2~=1 || sum(real(obj.chargeData.initial)~=...
				obj.chargeData.initial)>0 ||...
                sum(obj.chargeData.initial-obj.chargeData.minimum<=0)>0
                error(['chargeData.initial must have a real non-negative ',...
					'value for each receiver. Each value must be greater ',...
					'than chargeData.minimum.']);
            end
            
            [n1,n2] = size(obj.chargeData.threshold);
            if n1~=obj.nr || n2~=1 || sum(real(obj.chargeData.threshold)~=...
				obj.chargeData.threshold)>0 ||...
                sum(obj.chargeData.threshold-obj.chargeData.minimum<0)>0
                error(['chargeData.threshold must have a real non-negatve ',...
				'value for each receiver. Each value must be at least ',...
				'chargeData.minimum.']);
            end

            [n1,n2] = size(obj.chargeData.maximum);
            if n1~=obj.nr || n2~=1 || sum(real(obj.chargeData.maximum)~=...
				obj.chargeData.maximum)>0 ||...
                sum(obj.chargeData.maximum-obj.chargeData.threshold<0)>0
                error(['chargeData.maximum must have a real positive value ',...
					'for each receiver. Each value must be at least the ',...
					'corresponding chargeData.threshold']);
            end
            
            [n1,n2] = size(obj.constraints.maxCurr);
            if n1~=obj.n || n2~=1 || sum(real(obj.constraints.maxCurr)~=...
				obj.constraints.maxCurr)>0 || sum(obj.constraints.maxCurr<=0)>0
                error('maxCurr must have a real positive value for each device.');
            end

            %[n1,n2] = size(obj.constraints.maxPapp);
            %if n1~=1 || n2~=1 || real(obj.constraints.maxPapp)~=...
			%	obj.constraints.maxPapp || obj.constraints.maxPapp<=0
            %    error('maxPapp must be a real non-negative scalar.');
            %end
            
            [n1,n2] = size(obj.constraints.maxPact);
            if n1~=1 || n2~=1 || real(obj.constraints.maxPact)~=...
				obj.constraints.maxPact || obj.constraints.maxPact<=0
                error('maxPact must be a real non-negative scalar.');
            end
            
            [n1,n2] = size(obj.deviceData);
            if n1~=obj.nr || n2~=1
                error('Only one DeviceData object for each device is allowed');
            end
            for i=1:obj.nr
                if ~obj.deviceData(i).check()
                    error(['Invalid deviceData at index ',num2str(i)]);
                end
                [~,max_input_current] = obj.deviceData(i).domain_convACDC();
                [min_charge_current,~] = obj.deviceData(i).domain_effectiveChargeCurrent();
                if max_input_current < obj.constraints.maxCurr(obj.nt + i)
                    error(['The informed maxCurr for device ', num2str(i), ' is out of domain.']);
                end
                if min_charge_current > -maxId(i)
                    error(['The maximum requested discharge current for device ',...
                        num2str(i), ' is out of domain.']);
                end
            end
        end
        %Verify if a given solution is valid. Possible result values:
        %0: valid
        %1: incompatible size
        %2: unreached charge threshold
        %3: very high current
        %4: very high power (apparent)
        %5: very high power (active)
        %6: offline device (charge dropped below the minimum)
        %solution: ntxtime matrix with the transmitting voltages of each time slot
        function [result, QLog] = verify(obj, solution)
            obj = check(obj);

            %charge vector
            q = obj.chargeData.initial;
            QLog = q;

			[nt, time] = size(solution);
			
			if time==0
                %it is only possible if the devices start charged
                if mean(obj.chargeData.initial>=obj.chargeData.threshold)==1
                    result = 0;%ok, valid
                else
                    result = 2;%unreached
                end
				return;
            end
			
            if nt~=obj.nt
                result = 1;
                return;
            end

            %integrating...
            for t=1:time
                %calculating the load resistance of each receiving device
                Rl = [];
                for r = 1:obj.nr
                    Rl = [Rl; obj.deviceData(r).getRLfromSOC(q(r)/...
                        obj.chargeData.maximum(r))];
                end
				
                %calculating the phasor current vector
                current = (obj.timeLine(t).Z+diag([zeros(nt,1);Rl]))\...
                    [solution(:,t);zeros(obj.nr,1)];

                %verifying some constraints
                if sum(abs(current)>obj.constraints.maxCurr)>0
                    result = 3;
                    return;
                end
                %if abs(current(1:obj.nt)'*solution(:,t))>obj.constraints.maxPapp
                %    result = 4;
                %    return;
                %end
                if real(current(1:obj.nt)'*solution(:,t))>obj.constraints.maxPact
                    result = 5;
                    return;
                end
                %converting the currents
                chargeCurrent = [];
                for r = 1:obj.nr
                    %the input current less the discharge current
                    curr = obj.deviceData(r).convACDC(abs(current(obj.nt+r)))-...
                        obj.timeLine(t).Id(r);
                    %the charge/discharge current
                    chargeCurrent = [chargeCurrent;...
                        obj.deviceData(r).effectiveChargeCurrent(curr)];
                end

                %updating the charge vector
                q = min(q + obj.dt*chargeCurrent, obj.chargeData.maximum);
                QLog = [QLog,q];

                if sum(q<=obj.chargeData.minimum)>0
					disp('Minimo');
					q - obj.chargeData.minimum
                    result = 6;%there is an offline device
                    return;
                end
            end

            if sum(q<obj.chargeData.threshold)>0
				disp('Threshold');
				q - obj.chargeData.threshold
                result = 2;%could not complete all charges
            else
                result = 0;%valid solution
            end
        end

        %Plot the chart progression for each device
        %individual is true if the charts are meant to be separated for
        %each device
        function result = plot(obj, solution, individual)
            [result, QLog] = obj.verify(solution);

            figure; hold on;
            
            curves = [];
            labels = {};
            for r = 1:obj.nr
                %the charge curve for this device
                p = plot(QLog(r,:),obj.plot_colors(mod(r,length(obj.plot_colors))));

                if individual
                    curves(end+1) = p;
                    labels{end+1} = ['Device ', num2str(r)];
                end
                
                %the minimum charge (Dashed line)
                yline(obj.chargeData.minimum(r),...
                    ['--', obj.plot_colors(mod(r,length(obj.plot_colors)))]);

                %the maximum charge (Dash-dot line)
                yline(obj.chargeData.maximum(r),...
                    ['-.', obj.plot_colors(mod(r,length(obj.plot_colors)))]);

                %the minimum charge (dotted line)
                yline(obj.chargeData.threshold(r),...
                    [':', obj.plot_colors(mod(r,length(obj.plot_colors)))]);

                if ~individual
                    title(['Device ', num2str(r)]);
                    legend('charge','minumum','maximum','threshold');
                    %create a new figure for the next device
                    if r~=obj.nr
                        figure;hold on;
                    end
                end
            end
            if individual
                legend(curves, labels);
            end
        end
		
		%IT HAS A HUGE PROBELM!! IT DOES NOT CONSIDER THE CAPACITY OF THE BATTERIES.
		%The FeasibleFuture represented as a hyper-dot cloud may lead to precision lacks due to compression of the
		%intermediate states. Indeed, as the FeasibleFuture space is often too large, storing every value with 
		%great precision would lead to excessive memory usage. Therefore, this function was designed for improving the
		%precision of a solution.
		% - solution: vector of structures with the following fields
		%	- approximated voltage vector
		%	- charge vector
		function [success, new_solution, n_iterations_list] = recover_voltage_progression(obj, solution, max_iterations)
		
			success = true; %default
			n_iterations_list = []; %the number of iterations of each fine_adjustment call
			
			new_solution.Q = [];
			new_solution.V = [];
			
			[nt, n_slots] = size(solution.V);
			[nr, n_slots1] = size(solution.Q);
			
			if nt~=obj.nt || nr~=obj.nr || n_slots1~=n_slots
				error('Invalid sizes for the solution!!!');
			end
			
			%previous charge
			Q0 = obj.chargeData.initial;
			
			for i=1:n_slots
				
				Q = solution.Q(:,i);
				
				Ic = (Q - Q0)/obj.dt; %the required effective charge current
				
				target_reacheable = true; %default
				
				%some functions have constant intervals in their domain, so their inverse is not a function.
				%However, as they are monotonically increasing, we can manage the inverse as a function whose
				%image is an interval.
				min_targetIr = zeros(obj.nr,1);
				max_targetIr = zeros(obj.nr,1);
				
				RL = zeros(obj.nr, 1);
				
				for r = 1:obj.nr
					%inferring the load resistance for Q0
					SOC = Q0(r)/obj.chargeData.maximum(r);
					if SOC < 0 || SOC >1
						disp('Invalid SOC!');
						success = false;
						return;
					end
					RL(r) = obj.deviceData(r).getRLfromSOC(SOC);
				
					[min_Ic, max_Ic] = obj.deviceData(r).domain_iEffectiveChargeCurrent();
					%is this effective charge current possible?
					if min_Ic <= Ic(r) && Ic(r) <= max_Ic
						%the required input DC current
						[In0, In1] = obj.deviceData(r).iEffectiveChargeCurrent(Ic(r));
						
						In0 = In0 + obj.timeLine(i).Id(r);% - NPortPowerProblems.security_gap;
						In1 = In1 + obj.timeLine(i).Id(r);% + NPortPowerProblems.security_gap;
						
						[min_In, max_In] = obj.deviceData(r).domain_iConvACDC();
						
						if In1-In0 < NPortPowerProblems.tolerance
							%In1 = In0
							In = (In1 + In0)/2;
							if In < min_In || In > max_In
								disp('Invalid input current!');
								success = false;
								return;
							else
								if abs(obj.deviceData(r).effectiveChargeCurrent(In - obj.timeLine(i).Id(r)) - Ic(r))...
									> NPortPowerProblems.tolerance
									disp('Device Data imprecision!!!');
									success = false;
									return;
								end
								%the required amplitude for the receiving current
								[min_targetIr(r), max_targetIr(r)] = obj.deviceData(r).iConvACDC(In);
								
								%verifying if the maximum current constraint is satisfied
								if min_targetIr(r) + NPortPowerProblems.tolerance > obj.constraints.maxCurr(obj.nt+r)
									disp('Very high min target ir!');
									succes = false;
									return;
								else
									max_targetIr(r) = min(max_targetIr(r), obj.constraints.maxCurr(obj.nt+r) - NPortPowerProblems.tolerance);
								end
							end
						else
							%the interest region of the domain
							min_In = max(min_In, In0);
							max_In = min(max_In, In1);
							if min_In <= max_In
								
								if abs(obj.deviceData(r).effectiveChargeCurrent(min_In - obj.timeLine(i).Id(r)) - Ic(r))...
									> NPortPowerProblems.tolerance ...
									|| abs(obj.deviceData(r).effectiveChargeCurrent(max_In - obj.timeLine(i).Id(r)) - Ic(r))...
									> NPortPowerProblems.tolerance
									disp('Device Data imprecision!!!');
									success = false;
									return;
								end
							
								%the required amplitude for the receiving current
								[min_targetIr(r), ~] = obj.deviceData(r).iConvACDC(min_In);
								[~, max_targetIr(r)] = obj.deviceData(r).iConvACDC(max_In);
								
								%verifying if the maximum current constraint is satisfied
								if min_targetIr(r) + NPortPowerProblems.tolerance > obj.constraints.maxCurr(obj.nt+r)
									disp('Very high min target ir (2)!');
									succes = false;
									return;
								else
									max_targetIr(r) = min(max_targetIr(r), obj.constraints.maxCurr(obj.nt+r) - NPortPowerProblems.tolerance);
								end
								
							else
								disp('Invalid In interval!');
								success = false;
								return;
							end
						end
					else
						success = false;
						return;
					end
					
					%a simple test just to be sure targetIr interval is ok
					coeff = rand;
					if abs(Ic(r) - obj.deviceData(r).effectiveChargeCurrent(...
						obj.deviceData(r).convACDC((1-coeff)*min_targetIr(r)+coeff*max_targetIr(r))-obj.timeLine(i).Id(r)...
						)) > NPortPowerProblems.tolerance
						error('DeviceData failure!!');
					end
					
					if abs(Ic(r) - obj.deviceData(r).effectiveChargeCurrent(...
						obj.deviceData(r).convACDC(max_targetIr(r))-obj.timeLine(i).Id(r)...
						)) > NPortPowerProblems.tolerance
						error('DeviceData failure!!');
					end
					
					if abs(Ic(r) - obj.deviceData(r).effectiveChargeCurrent(...
						obj.deviceData(r).convACDC(min_targetIr(r))-obj.timeLine(i).Id(r)...
						)) > NPortPowerProblems.tolerance
						error('DeviceData failure!!');
					end
					
				end
				
				%the impedance matrix for this timeSlot
				Z = obj.timeLine(i).Z + diag([zeros(obj.nt,1);RL]);
				
				%the constraints
				It = obj.constraints.maxCurr(1:obj.nt) - 2*NPortPowerProblems.tolerance;
				P = obj.constraints.maxPact - 2*NPortPowerProblems.tolerance;
				
				ttl = 1000;
				V = solution.V(:,i);
				
				while true
					
					%getting the current for V
					I = Z\[V; zeros(obj.nr,1)];
					it = I(1:obj.nt);
					ir = I(obj.nt+1:end);
					
					%guarantee the constraints are respected
					if V.'*real(it) > P
						k = sqrt(P/(V.'*real(it)));
						V = k*V;
						it = k*it;
						ir = k*ir;
					end
					
					for t=1:obj.nt
						if abs(it(t)) > It(t)
							k = It(t)/abs(it(t));
							V = k*V;
							it = k*it;
							ir = k*ir;
						end
					end
					
					for r=1:obj.nr
						if abs(ir(r)) > obj.constraints.maxCurr(obj.nt+r) - NPortPowerProblems.tolerance
							k = (obj.constraints.maxCurr(obj.nt+r) - NPortPowerProblems.tolerance)/abs(ir(r));
							V = k*V;
							it = k*it;
							ir = k*ir;
						end
					end
					
					%now for targetIr itself we must choose the vector inside the interval [min_targetIr, max_target_Ir]
					%which is the closest to the former receiving voltage vector, that is, abs(ir)
					
					%first: which currents already are inside the target interval?
					inside = abs(ir) <= max_targetIr & abs(ir) >= min_targetIr;
					%which ones are under the interval?
					under = abs(ir) < min_targetIr;
					%what about over?
					over = abs(ir) > max_targetIr;
					
					targetIr = inside.*abs(ir) + under.*min_targetIr + over.*max_targetIr;
					
					[dv, n_iterations] = fine_adjustment(Z, V, it, ir, It, targetIr, P, NPortPowerProblems.tolerance, max_iterations);
					n_iterations_list = [n_iterations_list; n_iterations];
					
					if isempty(dv)
						V = solution.V(:,i) + normrnd(0,1,obj.nt,1);
					else
						%the new voltage vector
						V = V + dv;
						
						%getting the current for the new V
						I = Z\[V; zeros(obj.nr,1)];
						
						if max(abs(I)-obj.constraints.maxCurr)> NPortPowerProblems.tolerance
							error('The obtained I vector is out of limits.');
						end
						
						ir = I(obj.nt+1:end);
						
						%will the next state be as planned?
						Q_ok = true;
						for r=1:obj.nr
							%getting the new charge
							in = obj.deviceData(r).convACDC(abs(ir(r))) - obj.timeLine(i).Id(r);
							q = Q0(r) + obj.dt * obj.deviceData(r).effectiveChargeCurrent(in);
							%if q is too far from the previous solution, below the minimum or below
							%the threshold and it is the last slot
							if abs(Q(r)-q) > 100*obj.dt*NPortPowerProblems.tolerance ||...
								q <= obj.chargeData.minimum(r) ||...
								(q < obj.chargeData.threshold(r) && i==n_slots)
								Q_ok = false;
								break;
							else
								Q(r) = q;
							end
						end
						
						if Q_ok
							%the next state is as planned
							new_solution.Q = [new_solution.Q, Q];
							new_solution.V = [new_solution.V, V];
							break; %this slot is ready
						else
							%no, so search for other solution
							V = solution.V(:,i) + normrnd(0,1,obj.nt,1);
						end
					end
					
					if ttl <=0 
						success = false;
						disp('TTL is down!');
						return;
					else
						ttl = ttl - 1;
					end
				end
				
				%use the calculated charge instead of the charge of the input solution to reduce errors
				Q0 = Q;
			end
		end
		
        %Dummie implementation
        function [solveable, solution] = solve(obj)
            solveable = false;
            solution = [];
        end
    end
end
