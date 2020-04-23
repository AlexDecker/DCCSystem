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
            if nt~=obj.nt
                result = 1;
                return;
            end

            if time==0
                %it is only possible if the devices start charged
                if mean(obj.chargeData.initial>=obj.chargeData.threshold)==1
                    result = 0;%ok, valid
                else
                    result = 2;%unreached
                end
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
                    result = 6;%there is an offline device
                    return;
                end
            end

            if sum(q<obj.chargeData.threshold)>0
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
		
		%The FeasibleFuture represented as a hyper-dot cloud may lead to precision lacks due to compression of the
		%intermediate states. Indeed, as the FeasibleFuture space is often too large, storing every value with 
		%great precision would lead to excessive memory usage. Therefore, this function was designed for improving the
		%precision of a solution.
		% - solution: vector of structures with the following fields
		%	- previous charge vector
		%	- approximated voltage vector
		%	- charge vector
		function new_solution = recover_voltage_progression(obj, solution, max_iterations)
		
			new_solution = [];
			
			for i=1:length(solution)
			
				slot.Q0 = solution(i).Q0;
				slot.Q = solution(i).Q;
				
				Ic = (slot.Q - slot.Q0)/obj.dt; %the required effective charge current
				
				target_reacheable = true; %default
				%some functions have constant intervals in their domain, so their inverse is not a function.
				%However, as they are monotonically increasing, we can manage the inverse as a function whose
				%image is an interval.
				min_targetIr = zeros(obj.nr,1);
				max_targetIr = zeros(obj.nr,1);
				
				RL = zeros(obj.nr, 1);
				
				for r = 1:obj.nr
					%inferring the load resistance for Q0
					SOC = slot.Q0(r)/obj.chargeData.maximum(r);
					RL(r) = obj.deviceData(r).getRLfromSOC(SOC);
				
					[min_Ic, max_Ic] = obj.deviceData(r).domain_iEffectiveChargeCurrent();
					%is this effective charge current possible?
					if min_Ic <= Ic(r) && Ic(r) <= max_Ic
						%the required input DC current
						[In0,In1] = deviceData(r).iEffectiveChargeCurrent(Ic(r));
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
					%the impedance matrix for this timeSlot
					Z = obj.timeLine(i).Z + diag(RL);
					V = solution(i).V;
					%getting the current for V
					I = Z\[V; zeros(obj.nr,1)];
					it = I(1:obj.nt);
					ir = I(obj.nt+1:end);
					%the constraints
					It = obj.constraints.maxCurr(1:obj.nt);
					P = obj.constraints.maxPact;
					tolerance = 1e-10;
					
					dv = fine_adjustment(Z, solution(i).V, it, ir, It, targetIr, P, tolerance, max_iterations);
					slot.V = solution(i).V + dv;
					new_solution(end+1) = slot;
				else
					%not feasible, give up
					new_solution = [];
					return;
				end
			end
		end
		
        %Dummie implementation
        function [solveable, solution] = solve(obj)
            solveable = false;
            solution = [];
        end
    end
end
