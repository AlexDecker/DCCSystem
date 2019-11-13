%To better undertanding the expected behavior of each method, see FeasiblePast.m
%This class is functional but very inefficient for a large number of devices.
%In special, for a large number of receiving devices, the algorithms in this class
%are also very unaccurated (unless a very large bagSize is used)

%the set of charge vectors is represented as a hash.
classdef FPDummie < FeasiblePast
    properties
        bag %set of charge vectors (hash)
        nLevels
        maxCharge
        bagSize
    end
    methods
        function obj = FPDummie(hashSize,nLevels,maxCharge)
            if hashSize<1
                error('hashSize must be at least 1');
            end
            if nLevels<2
                error('nLevels must be at least 2');
            end
            [s1,s2] = size(maxCharge);
            if s1==0 || s2~=1
                error('Unexpected dimensions for maxCharge');
            end
            obj.bag = cell(hashSize,1);%create the hash as a cell Array
            obj.nLevels = nLevels;
            obj.maxCharge = maxCharge;
            obj.bagSize = nLevels^length(maxCharge);
        end

        function new = newFeasiblePast(obj, targetSet, timeSlot, dt,...
            minCharge, maxCurr, maxPapp, maxPact, rlCellArray, convCellArray)
            %as this function is called from NPortChargingProblem, the arguments
            %are expected to be valid

            %create a new empty object
            new = FPDummie(length(obj.bag),obj.nLevels,maxCharge);
            %populate the new object with acceptable charge vectors
            for i=0:obj.bagSize-1
                %generate the lambda ([0,1]^nR) coefficient vector
                lambda = zeros(length(maxCharge),1);
                for j=1:length(maxCharge)
                    lambda(j) = mod(floor(i/obj.nLevels^(j-1)),obj.nLevels)/(obj.nLevels-1);
                end
                %generate the actual charge vector from lambda
                candidate = lambda.*threshold + (1-lambda).*maxCharge;
                %evaluate if it is a feasible past vector
                [isFPV, voltages, next] = isFPVector(candidate, targetSet, timeSlot,...
                    dt, minCharge, maxCurr, maxPapp, maxPact, rlCellArray, convCellArray);
                if isFPV
                    q.charge = candidate;%the charge values
                    
                    q.voltages = voltages;%the voltages which are able to charge a set of
                    %sevices until their voltages be equal to 'next'

                    q.next = next;%a charge vector from the target set

                    %inserting the new element
                    target = insert(target,q);
                end
            end
        end
        %tests if a candidate charge vector is able to reach any charge vector from target
        %set at the given timeslot. If so, the voltages applyied for this time slot are
        %returned with the reached charge vector from the target set ('next')
        function [isFPV, voltages, next] = isFPVector(candidate, targetSet, timeSlot,...
            dt, minCharge, maxCurr, maxPapp, maxPact, rlCellArray, convCellArray)
            %the load resistance for this initial charge
            n = length(timeSlot.Z);nr = length(obj.maxCharge);nt=n-nr;
            Rl = zeros(n,1);
            for r=1:nr
                Rl(nt+r) = interp(rlCellArray{r}(:,1),rlCellArray{r}(:,2),...
                    candidate(r)/obj.maxCharge(r));
            end
            %building the inverse of the impedance matrix
            iZ = eye(length(timeSlot.Z))/(timeSlot.Z+diag(Rl));

        end
        function target = generateTarget(obj,threshold,maxCharge)
            [s1,s2] = size(threshold);
            if s1~=length(obj.maxCharge) || s2~=1
                error('Unexpected dimensions for threshold');
            end
            [s1,s2] = size(maxCharge);
            if s1~=length(obj.maxCharge) || s2~=1
                error('Unexpected dimensions for maxCharge');
            end
            %create a new empty object
            target = FPDummie(length(obj.bag),obj.nLevels,maxCharge);
            %populate the new object with acceptable final charge vectors
            for i=0:obj.bagSize-1
                %generate the lambda ([0,1]^nR) coefficient vector
                lambda = zeros(length(maxCharge),1);
                for j=1:length(maxCharge)
                    lambda(j) = mod(floor(i/obj.nLevels^(j-1)),obj.nLevels)/(obj.nLevels-1);
                end
                %generate the actual charge vector from lambda
                q.charge = lambda.*threshold + (1-lambda).*maxCharge;
                %the rest of the fields are dummies
                q.voltages = [];
                q.next = [];
                %inserting the new element
                target = insert(target,q);
            end
        end

        function q = search(obj,chargeVector)
            %get the element from the hash
            h = hashIndex(obj,chargeVector);
            %search in the colision set
            for i=1:length(obj.bag{h})
                %if the vectors are equal
                if sum(obj.bag{h}(i).charge~=chargeVector)==0
                    q = obj.bag{h}(i);
                    return;
                end
            end
            %if not found
            q = [];
        end
        %
        %hash methods-----------------------------------------------
        function h = hashIndex(obj,chargeVector)
            [s1,s2] = size(chargeVector);
            if s1~=length(obj.maxCharge) || s2~=1
                error('Unexpected dimensions for chargeVector');
            end
            h = ceil(mean(chargeVector./obj.maxCharge)*length(obj.bag));
        end

        function obj = insert(obj,element)
            h = hashIndex(obj,element.charge);
            obj.bag{h} = [obj.bag{h},element];
        end
    end
end
