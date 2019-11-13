%The user might create a class inheriting from this one in order to use different
%algorithms based on the framework defined by NPortChargingProblem class. In order
%to allow an instance of the NPortChargingProblem to operate with the desired 
%FeasiblePast child, create an empty object and inform it as feasiblePastModel when
%builing the instance of NPortChargingProblem.

%The feasible past consists of a potentially infinit set of charging vectors which
%may reach at least one charge vector from other given set in a given moment.
classdef FeasiblePast
    methods
        %creates an empty object to be used as model in NPortChargingProblem
        function obj = FeasiblePast()
        end
        
        %create a new feasible past, that means, a set of all charge vectors
        %which are able to reach at least one chargeVector from targetSet in a
        %given timeSlot, whithout disrespecting minCharge, maxCurr, maxPapp and
        %maxPact constraints. See NPortChargingProblem class
        function new = newFeasiblePast(obj, targetSet, timeSlot, dt,...
            minCharge, maxCurr, maxPapp, maxPact, rlCellArray, convCellArray)
            new = [];
        end
        
        %generate a target set, that is, the set of all valid charge vectors for
        %the charging to be complete
        function target = generateTarget(obj, threshold, maxCharge)
            target = [];
        end
        
        %search a given charge vector in the set, returning a structure containing
        %the following fields:
        %   * charge: the chargeVector itself
        %   * voltages: the active voltage vector to turn q into the target
        %   * next: the charge vector from the target set
        function q = search(obj, chargeVector)
            q = [];
        end
    end
end
