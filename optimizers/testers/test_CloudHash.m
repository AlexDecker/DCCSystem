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
TIME_FAIL = zeros(2,0);

%a list of vectors inside the hash
s=0;
MEMO = zeros(nr,0);

%for memory usage sake, v = q.^2 and q0 = sqrt(q), so they
%do not need to be explicitly stored
MAX = 1000;
while s<MAX
	%insert a new vector-------------------------------
	%generating a new discretized vector
	if rand<0.1 & s>=1
		%try to insert an already inserted vector
		d = MEMO(:,randi(s));
	else
		%try to generate a effectively new vector
		d = round(rand(nr,1)*(nSegments-1));
	end
	%0.1% tolerance: the limits for the new charge vetor
	q0 = minQ+(d+0.001).*sides;
	q1 = minQ+(d+0.999).*sides;
	%generating a random vector based on d
	r = rand(nr,1);
	q = q0.*r+(1-r).*q1;
	%prepare the data to be inserted
	data.q = q;
	data.v = q.^2;
	data.q0 = sqrt(q);
	%Has the vector already been inserted?
	i = find(mean(MEMO==d*ones(1,s))==1);
	if isempty(i)
		MEMO = [MEMO,d];%insert into the memory
		tic;
		[cloud,success] = insert(cloud,data);
		TIME_OK = [TIME_OK;toc];
		if ~success
			error('Insertion error');
		end
		s = s+1;
		%verify if size increased
		if cloud.size()~=s
			error('The size did not change correctly after the insertion');
		end
	else
		%try to insert and get an error
		tic;
		[cloud,success] = insert(cloud,data);
		TIME_FAIL = [TIME_FAIL,[s;toc]];
		if success
			error('Success for inserting rather than the expected failure');
		end
		%verify if size increased
		if cloud.size()~=s
			error('The size changed after an insertion failure');
		end
	end
	
	%log the initial size of the hash in Mega Bytes
	w = whos('cloud');
	MB = [MB;w.bytes/1000000];
	if mod(s,ceil(MAX/100))==0
		disp([num2str(100*s/MAX),'%']);
	end
end

figure;
plot(MB);
ylabel('MB');
xlabel('Numeber of elements');
title('Size of the hash');

figure;
hold on;
plot(TIME_OK);
plot(TIME_FAIL(1,:),TIME_FAIL(2,:));
ylabel('s');
xlabel('Numeber of elements');
title('Insertion time');
legend('OK','FAIL');

L=[];
for i=1:hashSize
	L = [L;cloud.len(i)];
end
figure;
plot(L);
xlabel('Entry');
title('Number of elements per hash entry');