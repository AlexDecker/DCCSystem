clear all;
%Test file for CloudHash.
nr = 10;
minQ = rand(nr,1);
maxQ = minQ + rand(nr,1);
nSegments = 100;
hashSize = 1000;
cloud = CloudHash(hashSize, nSegments, minQ, maxQ);

%calculating each side of the minimal unity of charge
sides = (maxQ-minQ)/nSegments;
%thus, a vector discretized as d correspond to the q vectors
%resulting from the convex linear combination of the following
%vectors (excluding q0 itself)
%q0 = minQ+d.*sides
%q1 = minQ+(d+1).*sides

%functions to be verified: read, size, search and insert


MB = [];
TIME_OK = [];
TIME_FAIL = [];

%a list of vectors inside the hash
s=0;
MEMO = zeros(nr,0);

%for memory usage sake, v = q.^2 and q0 = sqrt(q), so they
%do not need to be explicitly stored

while s<50
	%insert a new vector-------------------------------
	ttl = 10000;
	while true
		ttl = ttl-1;%do not seach to the infinity
		%generating a new discretized vector
		d = round(rand(nr,1)*(nSegments-1))
		%Has the vector already been inserted?
		i = find(sum(MEMO==d*ones(1,s)));
		if isempty(i)
			MEMO = [MEMO,d];%insert into the memory
			s = s+1;
			break;%stop searching
		end
		if ttl==0
			error('I give up');
		end
	end
	%0.1% tolerance: the limits for the new charge vetor
	q0 = minQ+(d+0.001).*sides;
	q1 = minQ+(d+0.999).*sides;
	%generating a random vector based on d
	r = rand(nr,1);
	q = q0.*r+(1-r).*q1;
	%insert the value
	data.q = q;
	data.v = q.^2;
	data.q0 = sqrt(q);
	tic;
	[cloud,success] = insert(cloud,data);
	TIME_OK = [TIME_OK;toc];
	if ~success
		error('Insertion error');
	end
	
	%verify if size increased------------------------
	if cloud.size()~=s
		error('The size did not change correctly after the insertion');
	end
	
	%insert a vector that will lead to failure
	d = MEMO(:,randi(s));
	%0.1% tolerance: the limits for the new charge vetor
	q0 = minQ+(d+0.001).*sides;
	q1 = minQ+(d+0.999).*sides;
	%generating a random vector based on d
	r = rand(nr,1);
	q = q0.*r+(1-r).*q1;
	%insert the value
	data.q = q;
	data.v = q.^2;
	data.q0 = sqrt(q);
	tic;
	[cloud,success] = insert(cloud,data);
	TIME_FAIL = [TIME_FAIL;toc];
	if success
		error('Success for inserting rather than the expected failure');
	end
	
	%verify if size increased------------------------
	if cloud.size()~=s
		error('The size changed after an insertion failure');
	end
	
	%log the initial size of the hash in Mega Bytes
	w = whos('cloud');
	MB = [MB;w.bytes/1000000];
	disp('Complete iteration: SUCCESS');
end

figure;
plot(MB);
ylabel('MB');
xlabel('Numeber of elements');
title('Size of the hash');

figure;
hold on;
plot(TIME_OK);
plot(TIME_FAIL);
ylabel('s');
xlabel('Numeber of elements');
title('Insertion time');
legend('OK','FAIL');