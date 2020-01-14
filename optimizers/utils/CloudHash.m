%This data structure stores a n-dimension point cloud as a hash.
%The space is divided into homogeneous small sub-spaces (hyper-
%cubes), where at most one point can be stored.

classdef CloudHash
    properties(GetAccess=public,SetAccess=private)
	
        hash %used for quickly finding elements
        LEN %number of elements currently stored
        Q0 %list of previous charge vectors
        V %list of voltages
        nSegments
        minQ
        maxQ
		
    end
	
    methods(Access=public)
	
        %hashSize: number of positions for the hash
        %nSegments: number of segments
        %minQ: minimum values for each q position (exclusive)
        %maxQ: maximum values for each q position (inclusive)
		%So, minQ<Q<=maxQ
        %maxSize: maximum number of elements for the hash
        function obj = CloudHash(hashSize, nSegments, minQ, maxQ, maxSize)
			if ~isscalar(nSegments) | nSegments<10 | nSegments>255
				error('nSegments must be a scalar between 10 and 255');
			end
            obj.nSegments = nSegments;
			if sum(minQ<0)>0 | sum(minQ>maxQ)
				error('0<=minQ<=maxQ');
			end
            obj.minQ = minQ;
            obj.maxQ = maxQ;
			if ~isscalar(hashSize) | hashSize<1
				error('Invalid hashSize');
			end

            %create an empty hash
            obj.hash = [];
            for i=1:hashSize
                e.D = zeros(length(minQ),5*ceil(maxSize/hashSize),'uint8');
                e.INDEX = zeros(1,5*ceil(maxSize/hashSize),'uint32');
                e.LEN = 0;
                obj.hash = [obj.hash;e];
            end

            %initializating the arrays where the data is effectively stored
            obj.LEN = 0;
            obj.Q0 = zeros(length(minQ),maxSize);
            obj.V = zeros(length(minQ),maxSize);
        end

        %insert a value in the hash. return true if it was effectively
        %inserted, false if there was already an entry in the hash
        %close enough to q (thus it was not inserted).
        %q:actual charge vector
        %v:voltage vector
        %q0:previous charge vector
        function [obj,ret] = insert(obj,q,v,q0)
            [found,h,~,d] = obj.search(q);
            if found
				%do nothing
                ret = false;
            else
				%insert the new row into the arrays
                obj.LEN = obj.LEN+1;
                obj.Q0(:,obj.LEN) = q0;
                obj.V(:,obj.LEN)  = v;

                i = obj.hash(h).LEN+1; %the new index in the h hash entry
                if i>length(obj.hash(h).INDEX)%if out of bounds
                    %duplicate the allocated size
                    obj.hash(h).D = [obj.hash(h).D, obj.hash(h).D];
                    obj.hash(h).INDEX = [obj.hash(h).INDEX, obj.hash(h).INDEX];
                end
                %effectively insert the new entry
                obj.hash(h).D(:,i) = d;
                obj.hash(h).INDEX(i) = length(obj.LEN);
                obj.hash(h).LEN = i;

                ret = true;
            end
        end

        %get the number of elements inside the h-th hash entry
        function n = len(obj,h)
            n = obj.hash(h).LEN;
        end
		
		%get the total number of elements inside the object
		function s = size(obj)
			s = 0;
			for h=1:length(obj.hash)
				s = s+obj.len(h);
			end
            if s~=obj.LEN
                error('Uncompatible data about the size of the hash');
            end
		end

        %read an entry given the indices hash(h) and inside the hash entry(i)
        function [D,Q0,V] = read(obj, h, i)
            D = obj.hash(h).D(:,i);
            Q0 = obj.Q0(:,obj.hash(h).INDEX(i));
            V = obj.V(:,obj.hash(h).INDEX(i));
        end

        %search a given vector q in the hash. return the hash index, the index
		%inside the hash entry and the discretized form of q
        function [found,h,i,d] = search(obj,q)
            d = obj.discretize(q);
            %find the entry in the hash
            h = obj.hashFunction(d);
            %find the element in the hash entry
            l = obj.hash(h).LEN;
            if l>0
                i = find(mean(obj.hash(h).D(:,1:l)==d*ones(1,l))==1);
                found = ~isempty(i);
            else
                i=[];
                found=false;
            end
        end
    end
	
    methods(Access=public)
        
        %returns the representation of charge vector q considering
        %the sub-spaces indexing. The returned vector has always
        %natural positive values
        function d = discretize(obj,q)
            d = ceil(obj.nSegments*(q-obj.minQ)./(obj.maxQ-obj.minQ))-1;
        end

        %deterministic hash function based on polynomial rolling hash function
        function h = hashFunction(obj,d)
            v = 0;
            for i=1:length(d)
                v = v*97+d(i);
            end
            h = mod(v,length(obj.hash))+1;
        end
    end
end
