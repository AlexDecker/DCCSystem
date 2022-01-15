% The point clouds are represented as multivariated gaussian distributions.
% Each one of them is computed from a given state from the previous time-slot.
% The point clouds from the "union" are already in the queue to be computed.
% The candidate will be evaluated and if it provides enough new states it will
% also be included in the queue.
% To evaluate it, we calculate the probability of ...

clc;
clear all;

side  = 100;
sampleSize = 1000;
f = 10;
sheetUnion = zeros(side);
sheetCandidate = zeros(side);

population = 10;
mi = cell(0);
sigma = cell(0);
for i = 1:population
	mi{end+1} = randi(side, 2, 1);
	sigma{end+1} = rand(2); sigma{end} = side/f * (sigma{end}*sigma{end}.' + 2*eye(2));
	sample = round(mvnrnd(mi{end},sigma{end},sampleSize));
	sample = max(1,min(side, sample));
	for j = 1:sampleSize
		sheetUnion(sample(j,1), sample(j,2)) = 1;
	end
end

mi_ = randi(side, 2, 1);
sigma_ = rand(2); sigma_ = side/f * (sigma_*sigma_.' + 2*eye(2));
sample = round(mvnrnd(mi_,sigma_,sampleSize));
sample = max(1,min(side, sample));
for j = 1:sampleSize
	sheetCandidate(sample(j,1), sample(j,2)) = 1;
end

intersection = sum(sum(sheetCandidate & sheetUnion));

candidate = sum(sum(sheetCandidate));

disp('Probability of a future state chosen through the candidate branch to be chosen through the other branches:');
if candidate == 0
	disp('100%'); %default
else
	disp([num2str(100*intersection/candidate),'%']);
end

image(20*(2*sheetCandidate + sheetUnion))