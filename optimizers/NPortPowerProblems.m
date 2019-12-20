%the version of the problems where the solver is aware about all parameters in the
%system. Minimize the charge time for a given WPT system or define a voltage progression
%for ensuring all devices ramain active for at least maxTau units of time
classdef NPortPowerProblems
    properties
        feasibleFutureModel %used to allow multiple algorithms to manage the feasibleFutures

        %time variant data
        timeLine
        dt
        
        %charge constraints and specific data
        chargeData

        %device behavior specific data
        deviceData

        %general constraints (not regarding charge)
        constraints
    end
    properties(Access=private)
        n %number of devices
        nt %number of transmitting devices
        nr %number of receiving devices
        maxTau %number of time slots
    end
    methods
        %timeLine: data for each time slot. Vector of structs with the following fields
        %   *Z: Impedance matrix without considering saturation or charge-dependent resistance
        %   *Id: discharge current (DC, real scalar)
        %dt: duration of one timeslot (real scalar)
        %chargeData: struct with vectors containing data about charge levels
        %   *minimum: vector with the minimum value of charge each device can handdle.
        %   *initial: vector with the initial charge of each device
        %   *threshold: vector with the minimum charge of each device at the end
        %   *maximum: the maximum charge supported by each battery
        %deviceData: DeviceData object list (one for each receiver)
        %   *rlCellArray: Lookup table for load resistance given the SOC
        %   *convCellArray: Lookup table for DC current given the AC input amplitude,
        %   from 0A to maxCurr
        %   *chargeCellArray: Lookup table for the effective charge current givent the
        %   dc input (from 0A to tha max converted current)
        %constraints: struct contaning
        %   *maxCurr: vector which contains the maximum current amplitude for each device
        %   *maxPapp: maximum allowed apparent power
        %   *maxPact: maximum allowed active power
        %feasibleFutureModel: an empty object from a class which inherits from FeasibleFuture which
        %   is used to create new objects and manage the different algorithms according to user's will
        function obj = NPortPowerProblems(timeLine, dt, chargeData, deviceData,...
            constraints, feasibleFutureModel)
            obj.timeLine = timeLine;
            obj.dt = dt;
            obj.chargeData = chargeData;
            obj.deviceData = deviceData;
            obj.constraints = constraints;
            obj.feasibleFutureModel = feasibleFutureModel;
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
                error('Id must be a column vector with one value for each receiving device');
            end
            obj.nr = n1;%receiving
            obj.nt = obj.n-obj.nr; %transmitting

            %verifying the timeline
            [n1,n2] = size(obj.timeLine);
            if n1==0 || n2~=1
                error('Invalid format for timeLine');
            end
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
                if sum(sum(real(obj.timeLine(t).Z-diag(diag(obj.timeLine(t).Z)))~=0))>0
                    error(['No real values outside the main diagonal. t=',num2str(t)]);
                end
                if sum(imag(diag(obj.timeLine(t).Z))~=0)>0
                    error(['Resonance is required. t=',num2str(t)]);
                end
                if sum(sum(obj.timeLine(t).Z-obj.timeLine(t).Z.'~=0))>0
                    error(['The mutual inductance must be symmetrical. t=',num2str(t)]);
                end
                if sum(abs(obj.timeLine(t).Id)~=obj.timeLine(t).Id)>0
                    error(['The discharge currents must be real positive. t=',num2str(t)]);
                end
            end
            
            [n1,n2] = size(obj.dt);
            if n1~=1 || n2~=1 || real(obj.dt)~=obj.dt || obj.dt<=0
                error('dt must be a real positive escalar.');
            end
            
            [n1,n2] = size(obj.chargeData.minimum);
            if n1~=obj.nr || n2~=1 || sum(real(obj.chargeData.minimum)~=obj.chargeData.minimum)>0 ||...
                sum(obj.chargeData.minimum<0)>0
                error('chargeData.minimum must have a real non-negative value for each receiver.');
            end

            [n1,n2] = size(obj.chargeData.initial);
            if n1~=obj.nr || n2~=1 || sum(real(obj.chargeData.initial)~=obj.chargeData.initial)>0 ||...
                sum(obj.chargeData.initial-obj.chargeData.minimum<=0)>0
                error(['chargeData.initial must have a real non-negative value for each receiver.',...
                    'Each value must be greater than the chargeData.minimum.']);
            end
            
            [n1,n2] = size(obj.chargeData.threshold);
            if n1~=obj.nr || n2~=1 || sum(real(obj.chargeData.threshold)~=obj.chargeData.threshold)>0 ||...
                sum(obj.chargeData.threshold-obj.chargeData.initial<0)>0
                error(['chargeData.threshold must have a real non-negatve value for each receiver. ',...
                    'Each value must be at least the chargeData.initial.']);
            end

            [n1,n2] = size(obj.chargeData.maximum);
            if n1~=obj.nr || n2~=1 || sum(real(obj.chargeData.maximum)~=obj.chargeData.maximum)>0 ||...
                sum(obj.chargeData.maximum-obj.chargeData.threshold<0)>0
                error(['chargeData.maximum must have a real positive value for each receiver. ',...
                    'Each value must be at least the corresponding chargeData.threshold']);
            end
            
            [n1,n2] = size(obj.constraints.maxCurr);
            if n1~=obj.n || n2~=1 || sum(real(obj.constraints.maxCurr)~=obj.constraints.maxCurr)>0 ||...
                sum(obj.constraints.maxCurr<=0)>0
                error('maxCurr must have a real positive value for each device.');
            end

            [n1,n2] = size(obj.constraints.maxPapp);
            if n1~=1 || n2~=1 || real(obj.constraints.maxPapp)~=obj.constraints.maxPapp ||...
                obj.constraints.maxPapp<=0
                error('maxPapp must be a real non-negative scalar.');
            end
            
            [n1,n2] = size(obj.constraints.maxPact);
            if n1~=1 || n2~=1 || real(obj.constraints.maxPact)~=obj.constraints.maxPact ||...
                obj.constraints.maxPact<=0
                error('maxPact must be a real non-negative scalar.');
            end
            
            [n1,n2] = size(obj.deviceData);
            if n1~=obj.nr || n2~=1
                error('There must be one DeviceData object for each receiving device');
            end
            for i=1:obj.nr
                if ~obj.deviceData(i).check()
                    error(['Invalid deviceData at index ',num2str(i)]);
                end
            end
        end
        %Verify if a given solution is valid
        %Result values:
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
            if nt~=obj.nt || time==0
                result = 1;
                return;
            end
            
            %integrating...
            for t=1:time
                %calculating the load resistance of each receiving device
                Rl = [];
                for r = 1:obj.nr
                    Rl = [Rl; obj.deviceData(r).getRLfromSOC(q(r,end)/obj.chargeData.maximum(r))];
                end       
                %calculating the phasor current vector
                current = (obj.timeLine(t).Z+diag([0;0;Rl]))\[solution(:,t);zeros(obj.nr,1)];
                %verifying some constraints
                if sum(abs(current)>obj.constraints.maxCurr)>0
                    result = 3;
                    return;
                end
                if abs(current(1:obj.nt)'*solution(:,t))>obj.constraints.maxPapp
                    result = 4;
                    return;
                end
                if real(current(1:obj.nt)'*solution(:,t))>obj.constraints.maxPact
                    result = 5;
                    return;
                end
                %converting the currents
                chargeCurrent = [];
                for r = 1:obj.nr
                    %the input current less the discharge current
                    curr = obj.deviceData(r).convACDC(abs(current(obj.nt+r)))-obj.timeLine(t).Id(r);
                    %the charge/discharge current
                    chargeCurrent = [chargeCurrent; obj.deviceData(r).effectiveChargeCurrent(curr)];
                end

                %updating the charge vector
                q = max(...
                        min(q + obj.dt*chargeCurrent, obj.chargeData.maximum),...
                        obj.chargeData.minimum);
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
        %Solves a given instance of the N-Port charging problem with up to maxTau time slots
        function [solveable, solution] = solveCharging(obj)
            %some verifications
            obj = check(obj);
            solution = [];

            %the initial state is already a valid solution?
            if mean(obj.chargeData.initial >= obj.chargeData.threshold)==1
                solveable = true;
                return;
            end
            
            %the initial state is invalid?
            if mean(obj.chargeData.initial < obj.chargeData.minimal)==1
                solveable = false;
                return;
            end

            %the final charge set (which is reached if the charging process is successful)
            initialSet = generateInitial(obj.feasibleFutureModel, obj.chargeData);

            %create the feasible futures up to the maximum number of time slots
            fFutureList = initialSet;
            for t=1:obj.maxTau 
                %this function creates a set with the states which are reacheable 
                [finalVector, fFuture] = newfeasibleFuture(obj.feasibleFutureModel,...
                    fFutureList(end),obj.timeLine(t), obj.dt, obj.chargeData,...
                    obj.deviceData, obj.constraints);
                if ~isempty(finalVector)
                    %solution found. building the voltage progression
                    solution = finalVector;
                    for i=length(fFutureList):-1:2 %the first element is the initial state
                        %search for the previous element
                        q = search(fFutureList(i),solution(1).previous);
                        %add a column
                        solution = [q.voltages, solution];
                    end
                    solveable = true;
                    return;
                end
                %verify if there is at least one feasible state for the current time slot
                if isEmpty(fFuture)
                    break;%No, so there is no solution.
                end
                fFutureList = [fFutureList; fFuture];
            end
            %no solution found
            solveable = false;
        end

        %Solves a given instance of the N-Port power sourcing problem with exactly maxTau time slots
        function [solveable, solution] = solvePowerSourcing(obj)
            %some verifications
            obj = check(obj);
            solution = [];

            %the initial state is already a valid solution?
            if mean(obj.chargeData.initial >= obj.chargeData.threshold)==1
                solveable = true;
                return;
            end
            
            %the initial state is invalid?
            if mean(obj.chargeData.initial < obj.chargeData.minimal)==1
                solveable = false;
                return;
            end

            %the final charge set (which is reached if the charging process is successful)
            initialSet = generateInitial(obj.feasibleFutureModel, obj.chargeData);

            %create the feasible futures up to the maximum number of time slots
            fFutureList = initialSet;
            for t=1:obj.maxTau 
                %this function creates a set with the states which are reacheable 
                [finalVector, fFuture] = newfeasibleFuture(obj.feasibleFutureModel,...
                    fFutureList(end),obj.timeLine(t), obj.dt, obj.chargeData,...
                    obj.deviceData, obj.constraints);    
                %verify if there is at least one feasible state for the current time slot
                if isEmpty(fFuture)
                    break;%No, so there is no solution.
                end
                fFutureList = [fFutureList; fFuture];
            end
            
            if ~isempty(finalVector)
                %solution found. building the voltage progression
                solution = finalVector;
                for i=length(fFutureList)-1:-1:2 %the first element is the initial state
                    %search for the previous element
                    q = search(fFutureList(i),solution(1).previous);
                    %add a column
                    solution = [q.voltages, solution];
                end
                solveable = true;
            else
                %no solution found
                solveable = false;
            end
        end
    end
end