%This class abstracts an instance of the N-Port Power Sourcing Problem, that is,
%the N-Port Power Problem which aims to find a time-domain voltage function which
%guarantees no device will be disconnected due low battery along a specified time
%interval. 
classdef NPortSourcingProblem < NPortPowerProblems
    properties
    end
    methods
        function obj = NPortSourcingProblem(timeLine, dt, chargeData,...
            deviceData, constraints, feasibleFutureModel)
            obj@NPortPowerProblems(timeLine, dt, chargeData, deviceData,...
                constraints, feasibleFutureModel)
        end

        function [solveable, solution] = solve(obj)
			%some verifications
            obj = check(obj);
            solution.Q = [];
			solution.V = [];

            %the initial state is invalid?
            if mean(obj.chargeData.initial > obj.chargeData.minimum)~=1
                solveable = false;
                return;
            end

            %the initial charge set (unitary)
            initialSet = generateInitialSet(obj.feasibleFutureModel, obj.chargeData);

            %create the feasible futures up to the maximum number of time slots
            fFutureList = initialSet;
            for t=1:obj.maxTau
                if t==obj.maxTau
                    stop_if_threshold_reached = true;
                else
                    stop_if_threshold_reached = false;
                end
                %this function creates a set with the states which are reacheable
                %from the previous one
                [finalElement, fFuture] = newFeasibleFuture(obj.feasibleFutureModel,...
                    fFutureList(end),...
					obj.timeLine(t), obj.dt, obj.chargeData,...
                    obj.deviceData, obj.constraints, stop_if_threshold_reached);    

                %verify if there is at least one feasible state for the current 
                %time slot
                if isEmpty(fFuture)
                    solveable = false;
                    return;%No, so there is no solution.
                end

                fFutureList = [fFutureList; fFuture];
            end
            
            if ~isempty(finalElement)
                %solution found. building the voltage progression matrix
                element = finalElement;
                solution.V = element.voltage;
				solution.Q = element.charge;

                for i=length(fFutureList)-1:-1:2 %the first element is the initial 
                %state. Go back through the time slots annotating the employied
                %voltages
                    %search for the previous element
                    element = search(fFutureList(i),element.previous);
                    %add a column
                    solution.V = [element.voltage, solution.V];
					solution.Q = [element.charge, solution.Q];
                end

                solveable = true;
            else
                %no solution found
                solveable = false;
            end 
        end
    end
end
