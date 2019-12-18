%See FeasibleFuture class for more information.
%This class implements the feasible future as a large set of elements, each one
%being equivalent to a reachable state given an initial set of states and a given
%time slot.
classdef FFDummie < FeasibleFuture
    properties
        minDistance %the set understends two very close elements as the same
    end
    methods
        function obj = FFuture(minDistance)
            obj.minDistance = minDistance;
        end
        
        function new = newFeasibleFuture(obj, initialSet, timeSlot, dt,...
            chargeData, deviceData, constraints)
            while ~isRepresentative(new)
                %get any element
                i = choose(new, initialSet);
                %generate a next solution
                s = generateNextState(new, i, timeSlot, dt, chargeData, constraints);
                new = insert(new,s);%insert if s does not belong to new
            end
            new = calculateRL(new, deviceData);%the load resistances
        end
        
        function initial = generateInitialSet(obj, chargeData)
            
        end
        
        %search a given charge vector in the set, returning a structure containing
        %the following fields:
        %   * charge: the chargeVector itself
        %   * voltages: the active voltage vector to turn previous into q
        %   * previous: the charge vector from the initial set
        function q = search(obj, chargeVector)
            q = [];
        end

        %returns true if there is no element in the set
        function b = isEmpty(obj)
            b = true;
        end
    end
end
