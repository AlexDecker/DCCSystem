%the version of the problem where the solver is aware about all parameters in the
%system. Minimize the charge time for a given WPT system or decide if it is possible
%to finish the charging in less than t1 units of time.
classdef NPortChargingProblem
    properties
        timeLine
        
        dt

        minCharge
        initialCharge
        chargeThreshold
        maxCharge

        maxCurr
        maxPapp
        maxPact

        rlCellArray
        convCellArray

        n %number of devices
        nt %number of transmitting devices
        nr %number of receiving devices
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
        %rlCellArray: cell array containing lookup tables for the load resistance given the SOC
        %convCellArray: cell array containing lookup tables for the conversion efficiency given
        %   the input amplitude, from 0A to maxCurr
        %maxCurr: vector which contains the maximum current amplitude for each device
        %maxPapp: maximum allowed apparent power
        %maxPact: maximum allowed active power
        function obj = NPortChargingProblem(timeLine, dt, chargeData, rlCellArray, convCellArray,...
            maxCurr, maxPapp, maxPact)
            obj.timeLine = timeLine;
            obj.dt = dt;
            obj.minCharge = chargeData.minimum;
            obj.initialCharge = chargeData.initial;
            obj.chargeThreshold = chargeData.threshold;
            obj.maxCharge = chargeData.maximum;
            obj.rlCellArray = rlCellArray;
            obj.convCellArray = convCellArray;
            obj.maxCurr = maxCurr;
            obj.maxPapp = maxPapp;
            obj.maxPact = maxPact;
            obj = check(obj);
        end
        %check if this object is acceptable
        function obj = check(obj)
            
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
            
            [n1,n2] = size(obj.minCharge);
            if n1~=obj.nr || n2~=1 || sum(real(obj.minCharge)~=obj.minCharge)>0 ||...
                sum(obj.minCharge<0)>0
                error('minCharge must have a real non-negative value for each receiver.');
            end

            [n1,n2] = size(obj.initialCharge);
            if n1~=obj.nr || n2~=1 || sum(real(obj.initialCharge)~=obj.initialCharge)>0 ||...
                sum(obj.initialCharge-obj.minCharge<=0)>0
                error(['initialCharge must have a real non-negative value for each receiver.',...
                    'Each value must be greater than the minCharge.']);
            end
            
            [n1,n2] = size(obj.chargeThreshold);
            if n1~=obj.nr || n2~=1 || sum(real(obj.chargeThreshold)~=obj.chargeThreshold)>0 ||...
                sum(obj.chargeThreshold-obj.initialCharge<0)>0
                error(['chargeThreshold must have a real non-negatve value for each receiver. ',...
                    'Each value must be at least the initialCharge.']);
            end

            [n1,n2] = size(obj.maxCharge);
            if n1~=obj.nr || n2~=1 || sum(real(obj.maxCharge)~=obj.maxCharge)>0 ||...
                sum(obj.maxCharge-obj.chargeThreshold<0)>0
                error(['maxCharge must have a real positive value for each receiver. ',...
                    'Each value must be at least the corresponding chargeThreshold']);
            end
            
            [n1,n2] = size(obj.maxCurr);
            if n1~=obj.n || n2~=1 || sum(real(obj.maxCurr)~=obj.maxCurr)>0 ||...
                sum(obj.maxCurr<=0)>0
                error('maxCurr must have a real positive value for each device.');
            end

            [n1,n2] = size(obj.maxPapp);
            if n1~=1 || n2~=1 || real(obj.maxPapp)~=obj.maxPapp || obj.maxPapp<=0
                error('maxPapp must be a real non-negative scalar.');
            end
            
            [n1,n2] = size(obj.maxPact);
            if n1~=1 || n2~=1 || real(obj.maxPact)~=obj.maxPact || obj.maxPact<=0
                error('maxPact must be a real non-negative scalar.');
            end

            if length(obj.rlCellArray)~=obj.nr
                error('rlCellArray must have one element for each receiving device.');
            end
            for i=1:obj.nr
                [n1,n2] = size(obj.rlCellArray{i});
                if n1<2 || n2~=2
                    error(['Invalid rlCellArray element: i=',num2str(i)]);
                end
                if obj.rlCellArray{i}(1,1)~=0 || obj.rlCellArray{i}(end,1)~=1 ||...
                    sum(obj.rlCellArray{i}(:,2)<0)>0
                    error(['Invalid data for rlCellArray element: i=',num2str(i)]);
                end
                for j=2:n1
                    if obj.rlCellArray{i}(j,1)<obj.rlCellArray{i}(j-1,1) ||...
                        obj.rlCellArray{i}(j,2)<=obj.rlCellArray{i}(j-1,2)
                        error(['Invalid data for rlCellArray element: i=',num2str(i)]);
                    end
                end
            end
            
            if length(obj.convCellArray)~=obj.nr
                error('convCellArray must have one element for each receiving device.');
            end
            for i=1:obj.nr
                [n1,n2] = size(obj.convCellArray{i});
                if n1<2 || n2~=2
                    error(['Invalid convCellArray element: i=',num2str(i)]);
                end
                if obj.convCellArray{i}(1,1)~=0 || obj.convCellArray{i}(end,1)~=obj.maxCurr(obj.nt+i) ||...
                    sum(obj.convCellArray{i}(:,2)<0)>0
                    error(['Invalid data for convCellArray element: i=',num2str(i)]);
                end
                for j=2:n1
                    if obj.convCellArray{i}(j,1)<obj.convCellArray{i}(j-1,1) ||...
                        obj.convCellArray{i}(j,2)<obj.convCellArray{i}(j-1,2)
                        error(['convCellArray must be monotonic: i=',num2str(i)]);
                    end
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
            q = obj.initialCharge;
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
                    Rl = [Rl; interp1(obj.rlCellArray{r}(:,1),obj.rlCellArray{r}(:,2),...
                            q(r,end)/obj.maxCharge(r))];
                end       
                %calculating the phasor current vector
                current = (obj.timeLine(t).Z+diag([0;0;Rl]))\[solution(:,t);zeros(obj.nr,1)];
                %verifying some constraints
                if sum(abs(current)>obj.maxCurr)>0
                    result = 3;
                    return;
                end
                if abs(current(1:obj.nt)'*solution(:,t))>obj.maxPapp
                    result = 4;
                    return;
                end
                if real(current(1:obj.nt)'*solution(:,t))>obj.maxPact
                    result = 5;
                    return;
                end
                %converting the currents
                chargeCurrent = [];
                for r = 1:obj.nr
                    curr = interp1(obj.convCellArray{r}(:,1),obj.convCellArray{r}(:,2),...
                        abs(current(obj.nt+r)))-obj.timeLine(t).Id(r);
                    chargeCurrent = [chargeCurrent; curr];
                end

                %updating the charge vector
                q = min(q + obj.dt*chargeCurrent, obj.maxCharge);
                QLog = [QLog,q];

                if sum(q<obj.minCharge)>0
                    result = 6;%there is an offline device
                    return;
                end
            end

            if sum(q<obj.chargeThreshold)>0
                result = 2;%could not complete all charges
            else
                result = 0;%valid solution
            end
        end
        %Decides if a given instance of the N-Port charging problem can be solved
        %for a given time limit t1. Returns a possible solution if it exists
        function [solveable, solution] = solveDecisionVersion(obj, t1)
            solveable = true;
            solution = [];
        end
    end
end
