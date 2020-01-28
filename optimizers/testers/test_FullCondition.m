%This script tests a method for determinign if a charge vector cloud is approximatly full
%The test is based on counting the number of consecutive insertion failures due to repeated
%input vectors

%creating an empty cloud
nr = 5;
nt = 6;
minQ = rand(nr,1);
maxQ = minQ + rand(nr,1);
nSegments = 25;
hashSize = 1000;
MAX = 1000000;
cloud = CloudHash(hashSize, nSegments, minQ, maxQ, MAX, nt);

%generating new random vectors
failures = 0;
successes = 0;
consecutive_failures = 0;
max_consecutive = 0;
tic;
while true
    r = 0.999*rand(nr,1)+0.001;
    q = maxQ.*r + (1-r).*minQ;
    d = cloud.discretize(q);
    [found,~,~,~] = cloud.search(d);
    if found
        consecutive_failures = consecutive_failures + 1;
        failures = failures + 1;
        max_consecutive = [max_consecutive;max(max_consecutive(end),consecutive_failures)];
    else
        consecutive_failures = 0;
        successes = successes + 1;
        cloud = insert(cloud,d,rand(nt,1),rand(nr,1));
        if mod(successes,1000)==0
            disp(successes);
            toc
            tic;
        end
    end
end
