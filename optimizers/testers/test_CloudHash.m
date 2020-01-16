%Test file for CloudHash.
nr = 5;
nt = 6;
minQ = rand(nr,1);
maxQ = minQ + rand(nr,1);
nSegments = 100;
hashSize = 49999;
MAX = 100000;
cloud = CloudHash(hashSize, nSegments, minQ, maxQ, MAX, nt);

%discretization tools

V = rand(nt,nr);%projection matrix to easily verify the values of the returned v from q
minV = zeros(nt,1);
maxV = (nSegments-1)*V*ones(nr,1);
volt_conv = CloudHash(1, 255, minV, maxV, 0, nt);%just for discretizing v

Q0 = rand(nr,nr);%projection matrix to easily verify the values of the returned q0 from q
minQ0 = zeros(nr,1);
maxQ0 = (nSegments-1)*Q0*ones(nr,1);
q0_conv = CloudHash(1, nSegments, minQ0, maxQ0, 0, nt);%just for discretizing q0

%functions to be verified: read, size, search and insert

MB = [];
TIME_OK = [];
TIME_FAIL = zeros(2,0);
TIME_DUM = zeros(2,0);

%a list of vectors inside the hash
s=0;
MEMO = zeros(nr,0);

%for memory usage sake, v = V*d and q0 = Q0*d, so they
%do not need to be explicitly stored
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
	%generating a random vector based on d
    [q1,q2] = cloud.dediscretize(d);
	r = rand(nr,1)*0.999;%guarantee r~=1
	q = q1.*r+(1-r).*q2;
	%prepare the data to be inserted
    v = volt_conv.discretize(V*d);
	d0 = q0_conv.discretize(Q0*d);
	%Has the vector already been inserted?
	tic;
	i = find(mean(MEMO==d*ones(1,s))==1);
	TIME_DUM = [TIME_DUM,[s;(2+nt/nr)*toc]];%log spent time
    %(the penality is because Q0 and V do not have storage cost for dummie)

	if isempty(i)
		tic;
		MEMO = [MEMO,d];%insert into the memory
		TIME_DUM(2,end) = TIME_DUM(2,end)+toc;
		tic;
        [found,~,~,~] = cloud.search(d);
		cloud = insert(cloud,d,v,d0);
        final = toc;
		TIME_OK = [TIME_OK;final];
		if found
			error('Insertion error');
		end
		s = s+1;
		%verify if size increased
		if cloud.size()~=s
			error('The size did not change correctly after the insertion');
		end

        %log the size of the hash in Mega Bytes
        w = whos('cloud');
        MB = [MB;w.bytes/1000000];
        
	else
		%try to insert and get an error
		tic;
		[found,~,~,~] = cloud.search(d);
		TIME_FAIL = [TIME_FAIL,[s;toc]];
		if ~found
			error('Success for inserting rather than the expected failure');
		end
		%verify if size increased
		if cloud.size()~=s
			error('The size changed after an insertion failure');
		end
	end
	if mod(s,ceil(MAX/1000))==0
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
plot(TIME_DUM(1,:),TIME_DUM(2,:));
ylabel('s');
xlabel('Numeber of elements');
title('Insertion time');
legend('search+insert','search','dummie_search');

L=[];
for i=1:hashSize
	L = [L;cloud.len(i)];
end
figure;
plot(L);
xlabel('Entry');
title('Number of elements per hash entry');
figure;
histogram(L);
xlabel('Number of elements per entry');
title('Collisions histogram');
