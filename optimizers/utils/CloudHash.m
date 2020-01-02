%This data structure stores a n-dimension point cloud as a hash.
%The space is divided into homogeneous small sub-spaces (hyper-
%cubes), where at most one point can be stored.

classdef CloudHash
    properties(GetAccess=public,SetAccess=private)
	
        hash
        nSegments
        minQ
        maxQ
        weightVector
		
    end
	
    methods(Access=public)
	
        %hashSize: number of positions for the hash
        %nSegments: number of segments
        %minQ: minimum values for each q position (exclusive)
        %maxQ: maximum values for each q position (inclusive)
		%So, minQ<Q<=maxQ
        function obj = CloudHash(hashSize, nSegments, minQ, maxQ)
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
			%generating a weight vector for random projection
            r = rand(1,length(minQ));
			obj.weightVector = 10000*(hashSize-1)/(sum(r)*(nSegments-1))*r;
            %create an empty hash
            for i=1:hashSize
                obj.hash{i}.D = zeros(length(minQ),0,'uint8');
                obj.hash{i}.DATA = [];
            end
        end

        %insert a value in the hash. return true if it was effectively
        %inserted, false if there was already an entry in the hash
        %close enough to q (thus it was not inserted).
        %data: struct with the fields
			%.q:actual charge vector
            %.v:voltage vector
            %.q0:previous charge vector
        function [obj,ret] = insert(obj,data)
            [found,h,~,d] = obj.search(data.q);
            if found
				%do nothing
                ret = false;
            else
				%insert the new entry
                obj.hash{h}.D = [obj.hash{h}.D, d];
                obj.hash{h}.DATA = [obj.hash{h}.DATA, data];
                ret = true;
            end
        end

        %get the number of elements inside the h-th hash entry
        function n = len(obj,h)
            n = length(obj.hash{h}.DATA);
        end
		
		%get the total number of elements inside the object
		function s = size(obj)
			s = 0;
			for h=1:length(obj.hash)
				s = s+obj.len(h);
			end
		end

        %read an entry given the indices hash(h) and inside the hash entry(i)
        function data = read(obj, h, i)
            data = obj.hash{h}.DATA(:,i);
        end

        %search a given vector q in the hash. return the hash index, the index
		%inside the hash entry and the discretized form of q
        function [found,h,i,d] = search(obj,q)
            d = obj.discretize(q);
            %find the entry in the hash
            h = obj.hashFunction(d);
            %find the element in the hash entry
            i = find(mean(obj.hash{h}.D==d*ones(1,obj.len(h)))==1);
            found = ~isempty(i);
        end
    end
	
    methods(Access=private)
        
        %returns the representation of charge vector q considering
        %the sub-spaces indexing. The returned vector has always
        %natural positive values
        function d = discretize(obj,q)
            d = ceil(obj.nSegments*(q-obj.minQ)./(obj.maxQ-obj.minQ))-1;
        end

        %deterministic hash function based on Random projection.
        function h = hashFunction(obj,d)
            h = mod(round(obj.weightVector*d),length(obj.hash))+1;
        end
    end
end
