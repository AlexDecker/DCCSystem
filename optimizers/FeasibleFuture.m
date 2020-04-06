%The user might create a class inheriting from this one in order to use different
%algorithms based on the framework defined by NPortPowerProblems class. In order
%to allow an instance of the NPortPowerProblems to operate with the desired 
%FeasibleFuture child, create an empty object and inform it as feasibleFutureModel when
%builing the instance of NPortPowerProblems.

%The feasible future consists of a potentially infinit set of charging vectors which
%are reacheable by at least one charge vector from other given set in a given moment.
classdef FeasibleFuture
    methods
        %creates an empty object to be used as model in NPortPowerProblems
        function obj = FeasibleFuture()
        end
        
        %create a new feasible future, which means, a set of all charge vectors
        %which are reachable by at least one chargeVector from previousSet in a
        %given timeSlot, whithout disrespecting minCharge, maxCurr and
        %maxPact constraints. See NPortPowerProblems class. If there is an element
        %in 'new' whose charge vector is more then the threshold, it must be returned
        %as 'final'. If there is more than one element, the function must choose one
        %of them according to user's will. The stop_if_threshold_reached boolean
        %argument controls if the search for new feasible future charge vectors must
        %be kept even if a final vector was found.
        function [final,new] = newFeasibleFuture(obj, initialSet, timeSlot, dt,...
            chargeData, deviceData, constraints, stop_if_threshold_reached)
            new = [];
            final = [];
        end
        
        %generate an initial set, that is, the unitary set with the initial state
        function initial = generateInitialSet(obj, chargeData)
            initial = [];
        end
        
        %search a given discretized charge vector in the set, returning a structure
        %containing the following fields:
        %   * voltages: the active voltage vector to turn previous into q
        %   * previous: the discretized charge vector from the initial set
        function q = search(obj, chargeVector)
            q = [];
        end

        %returns true if there is no element in the set
        function b = isEmpty(obj)
            b = true;
        end
    end
end
