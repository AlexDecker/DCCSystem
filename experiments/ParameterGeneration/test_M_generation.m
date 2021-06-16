% This script tests if the proposed algorithm for generating valid couplig matrices
% produces realistic matrix samplings

clear all;
close all;
clc;

rng('shuffle')

% sample size (increase this value if you want more instances, use 0 if you already have enough ones)
sample_size = 0;%100;

% number of coils
n_coils = 6;

% coil dimensions (we admit only homogeneous systems)
R2 = 0.15; % external radius
R1 = 0.05;% internal radius
N = 25; % number of turns
pts = 750; % coil resolution
wire_radius = 0.001;
prototype = SpiralPlanarCoil(R2,R1,N,wire_radius,pts);

% simulation area restricted to a 1 m origin-centered cube
maxV = 0.5;
maxR = pi/2;

% Generating a matrix sample using neumann formula
for n = 1:sample_size
	group_list = [];
	
	% creating the coils as planar-spirals
	for c = 1:n_coils
		coil = translateCoil(prototype,unifrnd(-maxV,maxV),...
                                unifrnd(-maxV,maxV),unifrnd(-maxV,maxV));
        group.coils.obj = rotateCoilX(rotateCoilY(coil,unifrnd(-maxR,maxR)),unifrnd(-maxR,maxR));
		group.R = -1;group.C = -1;
		group_list = [group_list; group];
		
	end
	
	
	% creating an environment (encapsulates the impedance matrix and calculates
	% the couplings. Using vacuum magnetic permeability and a dummie frequency
	environment = Environment(group_list,1e+6,pi*4e-7);
	
	if check(environment)
		% calculate the mutual inductance matriz with no previous information
		environment = environment.evalM(-ones(n_coils) + diag(ones(n_coils,1)));
		
		M = environment.M;
		
		save([num2str(round(rand*10000000)),'.mat'],'M');
		
		disp(['Instance number ',num2str(n),' saved successfully']);
	end
	
end

% loading all availiable neumann-based matrices in this directory
n_instances = 0;
CONTROL = zeros(0,(n_coils^2-n_coils)/2);
CONTROL_ALL = [];
mat = dir('*.mat');
for q = 1:length(mat)
	disp(['Loading ',mat(q).name, '...']);
	load(mat(q).name);
	m = [];
	% vectorizing the inductance matrix
	for i = 1:n_coils
		for j = i+1:n_coils
			m = [m, M(i,j)];
		end
	end
	CONTROL = [CONTROL;m];
	n_instances = n_instances + 1;
	CONTROL_ALL = [CONTROL_ALL; m.'];
end

if isempty(CONTROL)
	error('No .mat file in this directory');
end

% pearson correlation matrix
P = corr(CONTROL,'type','spearman');
figure;
image(P,'CDataMapping','scaled');
title('Control correlation');

% removing outliers
CONTROL_ALL_ = CONTROL_ALL(CONTROL_ALL<mean(CONTROL_ALL)+3*std(CONTROL_ALL));

% fitting in a gamma distribution
[~,params] = gamfit(CONTROL_ALL_);
shape = (params(1,1)+params(2,1))/2
scale = (params(1,2)+params(2,2))/2

% by symmetry, all variables are identically distributed
x = linspace(min(CONTROL_ALL_),max(CONTROL_ALL_),1000);
figure;
hold on;
histogram(CONTROL_ALL_,50,'normalization','pdf')
axis manual
plot(x, gampdf(x,shape,scale));
title('Control Distribution');

% Generating a matrix sample using the proposed method
EXPERIMENT = [];
EXPERIMENT_ALL = [];

n = 0; % number of EXPERIMENT samples so far

while n < n_instances
	
	% the number of devices must be n_coils to match the CONTROL sample
	result.nt = randi(n_coils - 3) + 1;
	result.nr = n_coils - result.nt;
	
	% arbitrary
	result.timeLine_size = randi(10);
	result.nSegments = 5+randi(15);
	result.sample_size = 2+randi(10);

	% arbitrary values for the lookup tables
	result.deviceData = [];
	for r = 1:result.nr
		%random lookup tables for load resistance, current conversion and charge conversion
		[rlTable,convTable,chargeTable] = randomLookupTables();
		%manager for the lookup tables
		result.deviceData = [result.deviceData; DeviceData(rlTable,convTable,chargeTable)];
	end

	success=false;

	disp('Generating a feasible N-Port Power Problem instance');
	while ~success
		[success, solution, result.chargeData, result.constraints, result.timeLine, result.dt] = generateFeasibleNPPPInstance(...
			result.deviceData, result.nt, result.nSegments, result.timeLine_size, result.sample_size);
	end
	disp('DONE.');
	
	% for each time_slot
	for t = 1:result.timeLine_size
		% extracting the inductance matrix from the impedance matrix
		M = -imag(result.timeLine(t).Z) / 1e6;
		m = [];
		% vectorizing the inductance matrix
		for i = 1:n_coils
			for j = i+1:n_coils
				m = [m, M(i,j)];
			end
		end
		if n < n_instances
			EXPERIMENT = [EXPERIMENT; m];
			EXPERIMENT_ALL = [EXPERIMENT_ALL; m.'];
			n = n + 1;
		else
			break;
		end
	end
end

% pearson correlation matrix
P = corr(EXPERIMENT,'type','spearman');
figure;
image(P,'CDataMapping','scaled');
title('Experiment correlation');

% fitting in a gamma distribution
%[~,params] = gamfit(EXPERIMENT_ALL);
%shape = (params(1,1)+params(2,1))/2
%scale = (params(1,2)+params(2,2))/2

% by symmetry, all variables are identically distributed
x = linspace(min(EXPERIMENT_ALL),max(EXPERIMENT_ALL),1000);
figure;
hold on;
histogram(EXPERIMENT_ALL,50,'normalization','pdf')
axis manual
%plot(x, gampdf(x,shape,scale));
title('Experiment Distribution');

% Testando se possuem a mesma distribuição (ks multivariada?)

% Se não for da mesma distribuição, que tal selecionar uma sub-amostra usando
% metropolis-hastings?
%}