%This data structure stores a n-dimension point cloud as a hash.
%The space is divided into homogeneous small sub-spaces (hyper-
%cubes), where at most one point can be stored.

%For memory sake, only the discretized version of the vectors is stored
%so, the vectors are stored using uint8, which allows levels from 0 to 255

classdef CloudHash
    properties(Constant)
        collision_factor = 3 %factor for the maximum allowed number of collisions
        initial_pool_size = 100 %initial size for the pools
    end
    properties(GetAccess=public,SetAccess=private)
        %sizes
        s %number of entries
        c %maximum number of collisions inside the hash
        nr %number of positions of the key
        nt %number of positions of V
        
        %hashes
        D %the hash with the key (n x s*c)
        LEN %a vector with the number of elements per hash entry
        D0 %has the same structure of hash, but stores the DO vectors
        V  %has the same structure of hash, but stores the V vectors
        A  %boolean value attached to each d vector

        %to handle escessive collisions
        pn %number of elements in the POOL
        ps %maximum number of elements in the POOL
        POOL_D %suppose index d
        POOL_D0 %the vector POOL_Q0(:,d) refers to d
        POOL_V %POOL_V(:,id) as well
        POOL_A %as well

        %for building the discretized version of q
        nSegments
        minQ
        maxQ
    end
	
    methods(Access=public)
	
        %nSegments: number of segments
        %minQ: minimum values for each q position (exclusive)
        %maxQ: maximum values for each q position (inclusive)
		%So, minQ<Q<=maxQ
        %nt: size of V
        function obj = CloudHash(hashSize, nSegments, minQ, maxQ, maxSize, nt)
            obj.nSegments = nSegments;
            obj.minQ = minQ;
            obj.maxQ = maxQ;

            obj.s = hashSize;
            obj.c = obj.collision_factor*ceil(maxSize/hashSize);
            obj.nr = length(minQ);
            obj.nt = nt;

            %create an empty hash
            obj.D = zeros(obj.nr,obj.s*obj.c,'uint8');
            obj.LEN = zeros(obj.s,1,'uint8');
            obj.D0 = zeros(obj.nr,obj.s*obj.c,'uint8');
            obj.V = zeros(obj.nt,obj.s*obj.c,'uint8');
            obj.A = false(obj.s,obj.c);

            %initializating the pools
            obj.pn = 0;
            obj.ps = obj.initial_pool_size;
            obj.POOL_D = zeros(obj.nr,obj.initial_pool_size,'uint8');
            obj.POOL_D0 = zeros(obj.nr,obj.initial_pool_size,'uint8');
            obj.POOL_V = zeros(obj.nt,obj.initial_pool_size,'uint8');
            obj.POOL_A = false(1, obj.initial_pool_size);
        end
        
        %insert a value which is not already inside the hash
        function obj = insert(obj,d,v,d0)
            h = obj.hashFunction(d);
            obj.LEN(h) = obj.LEN(h)+1;
            if obj.LEN>obj.c
				%do not fit in the hash
                if obj.pn==obj.ps
                    %resize the pool
                    obj.POOL_D = [obj.POOL_D,zeros(obj.nr,obj.ps,'uint8')];
                    obj.POOL_D0 = [obj.POOL_D0,zeros(obj.nr,obj.ps,'uint8')];
                    obj.POOL_V = [obj.POOL_V,zeros(obj.nt,obj.ps,'uint8')];
                    obj.POOL_A = [obj.POOL_A,false(1,obj.ps)];
                    obj.ps = obj.ps*2;
                end
                obj.pn = obj.pn+1;
                obj.POOL_D(:,obj.pn) = d;
                obj.POOL_D0(:,obj.pn) = d0;
                obj.POOL_V(:,obj.pn) = v;
            else
				%insert the new elements into the right entry
                i = (h-1)*obj.c+double(obj.LEN(h));
                obj.D(:,i) = d;
                obj.D0(:,i) = d0;
                obj.V(:,i) = v;
            end
        end

        %get the number of elements inside the h-th hash entry
        function n = len(obj,h)
            n = obj.LEN(h);
        end
		
		%get the total number of elements inside the object
		function s = size(obj)
			s = sum(obj.LEN);
		end

        %read an entry given the indices hash(h) and inside the hash entry(j)
        function [D,D0,V] = read(obj, h, j)
            i = (h-1)*obj.c+j;
            D = obj.D(:,i);
            D0 = obj.D0(:,i);
            V = obj.V(:,i);
        end

        %read from the collision pool
        function [D,D0,V] = readFromPool(obj,j)
            D = obj.POOL_D(:,i);
            D0 = obj.POOL_D0(:,i);
            V = obj.POOL_V(:,i);
        end

        %search a given discretized vector d in the hash. return the hash index,
        %the index j inside the hash entry (if it is in the entry) or the index
        %pj of the pool if it is in the pool
        function [found,h,j,pj] = search(obj,d)
            %find the entry in the hash
            h = obj.hashFunction(d);
            %find the element in the hash entry
            l = obj.LEN(h);
            if l>0
                %calculating the limits of the entry
                i0 = (h-1)*obj.c+1;
                i1 = i0+double(obj.LEN(h))-1;
                j = find(mean(obj.D(:,i0:i1)==d*ones(1,l))==1);
                if isempty(j)
                    if l>obj.c %the pool is guaranteed to not be empty
                        %it may be in the pool
                        pj = find(mean(obj.POOL_D(:,1:obj.pn)==d*ones(1,obj.pn))==1);
                        found = ~isempty(pj);
                    else
                        pj = [];
                        found = false;
                    end
                else
                    pj = [];
                    found = true;
                end
            else
                %the entry is empty
                j=[];
                pj = [];
                found=false;
            end
        end

        %annotate a discrete vector d with 'true'
        function obj = annotate(obj,d)
            [found,h,j,pj] = obj.search(d);
            if found
                if isempty(pj)
                    obj.A(h,j)=true;
                else
                    obj.POOL_A(pj) = true;
                end
            end
        end

        %create a smaller structure with only the signed vectors
        function small = compact(obj)
            %calculate the number of elements of 'small'
            n = sum(uint32(obj.A))+sum(uint32(POOL_A));
            %new size for the hash
            s = ceil(obj.s/n);
            %creating the new object
            small = CloudHash(s,obj.nSegments,obj.minQ,obj.maxQ,n,obj.nt);
            %add the elements in the hash
            for h=1:obj.s
                for i=1:obj.LEN(h)
                    if obj.A(h,i)
                       %TODO 
                    end
                end
            end
            %add the elements in the pool
            for i=1:obj.pn
                if POOL_A(i)
                    small = small.insert(POOL_D(:,i),POOL_V(:,i),POOL_D0(:,i));
                end
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

        %generate the vectors whose convex linear combination q = x*q0+(1-x)*q1
        %for x in [0,1[ is the generalization of the set of vectors that can be
        %discretized as d
        function [q0,q1] = dediscretize(obj,d)
            %calculating each side of the minimal unity of charge
            sides = (obj.maxQ-obj.minQ)/obj.nSegments;
            %thus, a vector discretized as d correspond to the q vectors                                   
            %resulting from the convex linear combination of the following
            %vectors (excluding q0 itself)
            q0 = obj.minQ+d.*sides;
            q1 = obj.minQ+(d+1).*sides;
        end

        %deterministic hash function based on polynomial rolling hash function
        function h = hashFunction(obj,d)
            v = 0;
            for i=1:length(d)
                v = v*97+d(i);
            end
            h = mod(v,obj.s)+1;
        end
    end
end
