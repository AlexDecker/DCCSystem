%To better undertanding the expected behavior of each method, see FeasiblePast.m
%This class is functional but very inefficient for a large number of devices.
%In special, for a large number of receiving devices, the algorithms in this class
%are also very unaccurated (unless a very large bagSize is used)

%the set of charge vectors is represented as a hash.
classdef FPDummie < FeasiblePast
    properties(GetAccess=public,SetAccess=private)
        bag %set of charge vectors (hash)
        nLevels %number of discrete levels for the charge values
        multipliers %used to avoid hash collisions. See hashIndex
        minStoredCharges %minimal charges already in the bag
        maxStoredCharges %maximal charges already in the bag
    end
    methods
        function obj = FPDummie(hashSize,nLevels,nr)
            if hashSize<1
                error('hashSize must be at least 1');
            end
            if nLevels<2
                error('nLevels must be at least 2');
            end
            if nr<1
                error('nr must be at least 1');
            end
        
            obj.bag = cell(hashSize,1);%create the hash as a cell Array
            obj.nLevels = nLevels;
            
            %search for nr prime numbers greater than the number of discrate
            %levels of each charge value
            num = 2*nLevels;%start by the first 100 numbers
            while true
                p = primes(num);
                obj.multipliers = p(p>nLevels);
                if length(obj.multipliers)>=nr
                    obj.multipliers = obj.multipliers(1:nr);
                    break;
                end
                num = 5*num;%not enough, get exponencially more numbers (thats ok,
                %since this function is not very common to be executed when compared
                %to others and nr is usually small when compared to nLevels)
            end

            %we still do not have stored vectors, so these are dummie
            obj.minStoredCharges = NaN*ones(nr,1);
            obj.maxStoredCharges = NaN*ones(nr,1);
        end

        function new = newFeasiblePast(obj, targetSet, timeSlot, dt,...
            chargeData, deviceData, constraints)
            %as this function is called from NPortChargingProblem, the arguments
            %are expected to be valid
            
            %getting the number of receivers
            nr = length(chargeData.minimum);

            %create a new empty object
            new = FPDummie(length(obj.bag),obj.nLevels,nr);

            %All feasible past vectors pareto-dominates minLimit.
            %In this case, the maximum 
            minLimit = max(target.minStoredCharges - ...
                dt*constraints.maxCurrent(end-nr+1:end) + 
                dt*timeSlot.Id,...
                chargeData.minimum);
            %maxLimit pareto-dominates all feasible past vectors
            %Roughly, there is no charge current and still the largest values from
            %the target are reached
            maxLimit = min(target.maxStoredCharges + dt*timeSlot.Id, chargeData.maximum);

            %populate the new object with acceptable charge vectors
            bagSize = obj.nLevels^nr;
            for i=0:bagSize-1
                %generate the lambda ([0,1]^nR) coefficient vector
                lambda = zeros(nr,1);
                for j=1:nr
                    lambda(j) = mod(floor(i/obj.nLevels^(j-1)),obj.nLevels)/(obj.nLevels-1);
                end
                %generate the actual charge vector from lambda
                candidate = lambda.*minLimit + (1-lambda).*maxLimit;
                %evaluate if it is a feasible past vector
                [isFPV, voltages, next] = isFPVector(candidate, targetSet, timeSlot,...
                    dt, chargeData, deviceData, constraints);
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
            dt, chargeData, deviceData, constraints)
            %the load resistance for this initial charge
            n = length(timeSlot.Z);nr = length(chargeData.maxCharge);nt=n-nr;
            Rl = zeros(n,1);
            for r=1:nr
                Rl(nt+r) = interp(deviceData.rlCellArray{r}(:,1),deviceData.rlCellArray{r}(:,2),...
                    candidate(r)/chargeData.maxCharge(r));
            end
            %building the inverse of the impedance matrix
            iZ = eye(length(timeSlot.Z))/(timeSlot.Z+diag(Rl));

        end
        function target = generateTarget(obj,chargeData)
            %create a new empty object
            target = FPDummie(length(obj.bag),obj.nLevels,nr);
            
            %getting the number of receivers
            nr = length(chargeData.minimum);

            %populate the new object with acceptable charge vectors
            bagSize = obj.nLevels^nr;

            %populate the new object with acceptable final charge vectors
            for i=0:bagSize-1
                %generate the lambda ([0,1]^nR) coefficient vector
                lambda = zeros(length(chargeData.maximum),1);
                for j=1:length(chargeData.maximum)
                    lambda(j) = mod(floor(i/obj.nLevels^(j-1)),obj.nLevels)/(obj.nLevels-1);
                end
                %generate the actual charge vector from lambda
                q.charge = lambda.*chargeData.threshold + (1-lambda).*chargeData.maximum;
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
    end
    methods(Access=private)
        %hash methods-----------------------------------------------
        function h = hashIndex(obj,chargeVector)
            [s1,s2] = size(chargeVector);
            if s1~=length(obj.maxLimit) || s2~=1
                error('Unexpected dimensions for chargeVector');
            end
            %minLimit is 0, maxLimit is 1.
            normalized = (chargeVector-obj.minLimit)./(obj.maxLimit-obj.minLimit);
            %transforms into a value between 1 and hashSize (eventually getting anormal values)
            h = ceil(mean(normalized)*length(obj.bag));
            %verify eventual anormal values
            if h<=0
                h=1;
            else
                if h>length(obj.bag)
                    h=length(obj.bag);
                end
            end
        end

        function obj = insert(obj,element)
            h = hashIndex(obj,element.charge);
            obj.bag{h} = [obj.bag{h},element];
            %update the limits for the stored vectors
            obj.minStoredCharges = min(obj.minStoredCharges,element.charge);
            obj.maxStoredCharges = max(obj.maxStoredCharges,element.charge);
        end
    end
end
