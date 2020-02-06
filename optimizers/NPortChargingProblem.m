%This class abstracts an instance of the N-Port Charging Problem, that is,
%the N-Port Power Problem which aims to minimize the time required for all
%devices to exceed a given threshold.
classdef NPortChargingProblem < NPortPowerProblems
    properties
    end
    methods
        function obj = NPortChargingProblem(timeLine, dt, chargeData,...
            deviceData, constraints, feasibleFutureModel)
            obj@NPortPowerProblems(timeLine, dt, chargeData, deviceData,...
                constraints, feasibleFutureModel)
        end

        function [solveable, solution] = solve(obj)

            %some verifications
            obj = check(obj);

            solution = [];

            %the initial state is already a valid solution?
            if mean(obj.chargeData.initial >= obj.chargeData.threshold)==1
                solveable = true;
                return;
            end
            
            %the initial state is invalid?
            if mean(obj.chargeData.initial <= obj.chargeData.minimal)==1
                solveable = false;
                return;
            end

            %the initial charge set (unitary set)
            initialSet = generateInitialSet(obj.feasibleFutureModel, obj.chargeData);

            %create the feasible futures up to the maximum number of time slots
            fFutureList = initialSet;

            for t=1:obj.maxTau 
                %this function creates a set with the states which are reacheable 
                [finalElement, fFuture]=newfeasibleFuture(obj.feasibleFutureModel,...
                    fFutureList(end),obj.timeLine(t), obj.dt, obj.chargeData,...
                    obj.deviceData, obj.constraints);

                if ~isempty(finalElement)
                    %solution found. building the voltage progression
                    element = finalElement;
                    solution = element.voltage;

                    for i=length(fFutureList):-1:2 %the first element is the 
                    %initial state

                        %search for the previous element
                        element = search(fFutureList(i),element.previous);

                        %add a column
                        solution = [element.voltages, solution];
                    end

                    solveable = true;
                    return;
                end

                %verify if there is at least one feasible state for the current 
                %time slot
                if isEmpty(fFuture)
                    solveable = false;
                    return;%No, so there is no solution.
                end

                fFutureList = [fFutureList; fFuture];
            end

            %no solution found
            solveable = false;
        end
    end
end
