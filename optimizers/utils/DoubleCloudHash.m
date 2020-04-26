%This data structure stores a n-dimension point cloud as a hash.
%The space is divided into homogeneous small sub-spaces (hyper-
%cubes), where at most one point can be stored.

classdef DoubleCloudHash
    properties(Constant)
        collision_factor = 1 %factor for the maximum allowed number of collisions
        initial_pool_size = 100 %initial size for the pools
    end
    properties(GetAccess=public,SetAccess=private)
        %sizes
        s %number of entries
        c %maximum number of collisions inside the hash
        nr %number of positions of the key
        nt %number of positions of V
		maxSize %maximum number of elements

        multiplier %for the hash function
        
        %hashes
        D %the hash with the key (n x s*c)
		Q %the hash with the charge vectors
        LEN %a vector with the number of elements per hash entry
        D0 %has the same structure of hash, but stores the DO vectors
        V  %has the same structure of hash, but stores the V vectors
        A  %boolean value attached to each d vector

        %to handle excessive collisions
        pn %number of elements in the POOL
        ps %maximum number of elements in the POOL
        POOL_D %suppose index d
		POOL_Q %the vector POOL_Q(:,d) refers to d
        POOL_D0 %the vector POOL_D0(:,d) refers to d
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
        function obj = DoubleCloudHash(hashSize, nSegments, minQ, maxQ, maxSize, nt)
            if nSegments<2
                error('DoubleCloudHash: the number of segments must be at least 2');
            end
            obj.nSegments = nSegments;
            obj.minQ = minQ;
            obj.maxQ = maxQ;
            
            %the actual hashSize is the largest prime less than or equal to hashSize
            %for better performance
            p = primes(hashSize);
            if length(p)>0
                hashSize = p(end);
            end
            %the multiplier of the hash function must also be the closest prime to
            %the number of segments
            p = primes(obj.nSegments);
            obj.multiplier = p(end);

            obj.s = hashSize;
            obj.c = ceil(obj.collision_factor*maxSize/hashSize);
            obj.nr = length(minQ);
            obj.nt = nt;
			obj.maxSize = maxSize;

            %create an empty hash
            obj.D = zeros(obj.nr,obj.s*obj.c,'uint8');
			obj.Q = zeros(obj.nr,obj.s*obj.c);
            obj.LEN = zeros(obj.s,1);
            obj.D0 = zeros(obj.nr,obj.s*obj.c,'uint8');
            obj.V = zeros(obj.nt,obj.s*obj.c);
            obj.A = false(obj.s,obj.c);

            %initializating the pools
            obj.pn = 0;
            obj.ps = obj.initial_pool_size;
            obj.POOL_D = zeros(obj.nr,obj.initial_pool_size,'uint8');
			obj.POOL_Q = zeros(obj.nr,obj.initial_pool_size);
            obj.POOL_D0 = zeros(obj.nr,obj.initial_pool_size,'uint8');
            obj.POOL_V = zeros(obj.nt,obj.initial_pool_size);
            obj.POOL_A = false(1, obj.initial_pool_size);
        end
        
        %insert a value which is not already inside the hash
        function obj = insert(obj,q,v,d0)
			d = obj.discretize(q);
            h = obj.hashFunction(d);
            obj.LEN(h) = obj.LEN(h)+1;
            if obj.LEN(h)>obj.c
				%do not fit in the hash
                if obj.pn==obj.ps
                    %resize the pool
                    obj.POOL_D = [obj.POOL_D,zeros(obj.nr,obj.ps,'uint8')];
					obj.POOL_Q = [obj.POOL_Q,zeros(obj.nr,obj.ps)];
                    obj.POOL_D0 = [obj.POOL_D0,zeros(obj.nr,obj.ps,'uint8')];
                    obj.POOL_V = [obj.POOL_V,zeros(obj.nt,obj.ps)];
                    obj.POOL_A = [obj.POOL_A,false(1,obj.ps)];
                    obj.ps = obj.ps*2;
                end
                obj.pn = obj.pn+1;
                obj.POOL_D(:,obj.pn) = uint8(d);
				obj.POOL_Q(:,obj.pn) = q;
                obj.POOL_D0(:,obj.pn) = uint8(d0);
                obj.POOL_V(:,obj.pn) = v;
            else
				%insert the new elements into the right entry
                i = (h-1)*obj.c+double(obj.LEN(h));
                obj.D(:,i) = uint8(d);
				obj.Q(:,i) = q;
                obj.D0(:,i) = uint8(d0);
                obj.V(:,i) = v;
            end
        end
		
		%get the total number of elements inside the object
		function s = countElements(obj)
			s = sum(obj.LEN);
		end
		
		%return true if the maximum capacity of the structure is reached
		function b = isFull(obj)
			b = obj.countElements() == obj.maxSize;
		end
		
        %read an entry given the indices hash(h) and inside the hash entry(j)
        function [D,Q,D0,V] = read(obj, h, j)
            i = (h-1)*obj.c+j;
            D = double(obj.D(:,i));
			Q = obj.Q(:,i);
            D0 = double(obj.D0(:,i));
            V = obj.V(:,i);
        end

        %read from the collision pool
        function [D,Q,D0,V] = readFromPool(obj,j)
            D = double(obj.POOL_D(:,j));
			Q = obj.POOL_Q(:,j);
            D0 = double(obj.POOL_D0(:,j));
            V = obj.POOL_V(:,j);
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
                %valid length of the entry
                len = min(double(l),obj.c);
                %calculating the limits of the entry
                i0 = (h-1)*obj.c + 1;
                i1 = i0 + len - 1;
                %searching in the hash entry
                j = find(mean(obj.D(:,i0:i1)==d*ones(1,len))==1);
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
            n = sum(uint32(obj.A))+sum(uint32(obj.POOL_A));
            %new size for the hash
            s = ceil(obj.s/n);
            %creating the new object
            small = DoubleCloudHash(s,obj.nSegments,obj.minQ,obj.maxQ,n,obj.nt);
            %add the elements in the hash
            for h=1:obj.s
                for i=1:obj.LEN(h)
                    if obj.A(h,i)
                       small = small.insert(obj.Q(:,i),obj.V(:,i),obj.D0(:,i));
                    end
                end
            end
            %add the elements in the pool
            for i=1:obj.pn
                if obj.POOL_A(i)
                    small = small.insert(obj.POOL_Q(:,i),obj.POOL_V(:,i),obj.POOL_D0(:,i));
                end
            end
        end

        %choose any charge vector to return
        function [Q,D] = any(obj)
            
            if obj.countElements()==0
                error('DoubleCloudHash.any: the object is empty')
            end

            h = randi(obj.s);%get any entry

            %get the next one until reaching a non-empty entry
            while obj.LEN(h)==0
                h = mod( h, obj.s ) + 1;
            end

            j = randi(obj.LEN(h));%get any element from the entry

            if j>obj.c
                %get any element from the pool
                [D,Q,~,~] = obj.readFromPool(randi(obj.pn));
            else
                %get the element from the hash
                [D,Q,~,~] = obj.read(h,j);
            end
        end
    end
	
    methods(Access=public)
        
        %returns the representation of charge vector q considering
        %the sub-spaces indexing. The returned vector has always
        %natural positive values, unless q is less than or equal
		%to minQ
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
                v = v*obj.multiplier+d(i);
            end
            h = mod(v,obj.s)+1;
        end
    end
end
