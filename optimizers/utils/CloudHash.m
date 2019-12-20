%This data structure stores a n-dimension point cloud as a hash.
%The space is divided into homogeneous small sub-spaces (hyper-
%cubes), where at most one point can be stored.

classdef CloudHash
    properties(GetAccess=public,SetAccess=private)
        hash
        nSegments
        minQ
        maxQ
        primeVector
    end
    methods(Static)
        
    end
    methods(Access=public)
        %hashSize: number of positions for the hash
        %nSegments: number of segments
        %minQ: minimum values for each q position
        %maxQ: maximum values for each q position
        function obj = CloudHash(hashSize, nSegments, minQ, maxQ)
            if ~isprime(hashSize)
                error('hash size must be prime');
            end
            obj.nSegments = nSegments;
            obj.minQ = minQ;
            obj.maxQ = maxQ;
            [obj,next] = generatePrimeVector(obj);
            if next>hashSize
                disp('Warning: hash size altered for consitency');
                hashSize = next;
            end
            %create an empty hash
            for i=1:hashSize
                obj.hash{i}.Q = [];
                obj.hash{i}.DATA = [];
            end
        end

        %insert a value in the hash. return true if it was effectively
        %inserted, false if there was already an entry in the hash
        %close enough to q (thus it was not inserted).
        %data: struct with the fields
            %.v:voltage vector
            %.q0:previous charge vector
        function [obj,ret] = insert(obj,q,data)
            [found,~,~] = obj.search(q);
            if ~found
                %insert the new entry
                obj.hash{h}.Q = [obj.hash{h}.Q, q];
                obj.hash{h}.DATA = [obj.hash{h}.DATA, data];
                ret = true;
            else
                %do nothing
                ret = false;
            end
        end

        %get the number of elements inside the h-th hash entry
        function n = len(obj,h)
            n = length(obj.hash{h}.Q);
        end

        %read an entry given the indices hash(h) and inside the hash entry(i)
        function [q,data] = read(obj, h, i)
            q = obj.hash{h}.Q(:,i);
            data = obj.hash{h}.DATA(:,i);
        end

        %search a given vector q in the hash
        function [found,q,data] = search(obj,q)
            d = obj.discretize(q);
            %find the entry in the hash
            h = obj.hashFunction(d);
            %find the element in the hash entry
            i = find(sum(obj.hash{h}.Q==q*ones(1,length(obj.hash{h}.Q))));
            if isempty(i)
                %return the found entry
                found = true;
                q = obj.hash{h}.Q(:,i);
                data = obj.hash{h}.DATA(:,i);
            else
                %not found
                found = false;
                q=[];data=[];
            end
        end
    end
    methods(Access=private)
        
        %returns the representation of charge vector q considering
        %the sub-spaces indexing. The returned vector has always
        %natural positive values
        function d = discretize(obj,q)
            d = floor(obj.nSegments*(q-obj.minQ)./(obj.maxQ-obj.minQ));
        end
        
        %generate a vector with the first length(q) prime numbers.
        %also generates a candidate for substituting hashSize if
        %needed.
        function [obj,next] = generatePrimeVector(obj)
            num = 100;
            while true
                p = primes(num);
                if length(p)>=length(obj.minQ)+1
                    break;
                else
                    num = num*10;
                end
            end
            obj.primeVector = p(1:length(obj.minQ));
            next = p(length(obj.minQ)+1);
        end

        %deterministic hash function based on Gobel numbering.
        %the input is a positive natural number vector.
        function h = hashFunction(obj,d)
            %the Hash Function:
            %residue of the dividion between sum(p_i^d_i) and hashSize
            %the p vector contains only unique prime entries.
            %the sum is unique for each vector d.
            h = 0;
            for i=1:length(d)
                %the residue of p_i^d_i/hashSize is equal to the residue
                %of p_i^r_i/hashSize, where r_i is the residue of 
                %d_i/(hashSize-1). This follows directly from Fermat's
                %little theorem
                r = mod(d(i),length(obj.hash)-1);
                %this function calculated the remainder
            end
        end
    end
end
