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
            solution = [];

            %the initial state is invalid?
            if mean(obj.chargeData.initial > obj.chargeData.minimal)~=1
                solveable = false;
                return;
            end

            %the initial charge set (unitary)
            initialSet = generateInitial(obj.feasibleFutureModel, obj.chargeData);

            %create the feasible futures up to the maximum number of time slots
            fFutureList = initialSet;
            for t=1:obj.maxTau 
                %this function creates a set with the states which are reacheable
                %from the previous one
                [finalElement, fFuture]=newfeasibleFuture(obj.feasibleFutureModel,...
                    fFutureList(end),obj.timeLine(t), obj.dt, obj.chargeData,...
                    obj.deviceData, obj.constraints);    

                %verify if there is at least one feasible state for the current 
                %time slot
                if isEmpty(fFiuture)
                    solveable = false;
                    return;%No, so there is no solution.
                end

                fFutureList = [fFutureList; fFuture];
            end
            
            if ~isempty(finalElement)
                %solution found. building the voltage progression matrix
                element = finalElement;
                solution = element.voltage;

                for i=length(fFutureList)-1:-1:2 %the first element is the initial 
                %state. Go back through the time slots annotating the employied
                %voltages
                    %search for the previous element
                    element = search(fFutureList(i),element.previous);
                    %add a column
                    solution = [element.voltage, solution];
                end

                solveable = true;
            else
                %no solution found
                solveable = false;
            end 
        end
    end
end
