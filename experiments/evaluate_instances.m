clc;
clear all;

n_files = 100;
folder = 'power_sourcing';
prefix = 'exp_';

found_solutions = 0;
failures = 0;
errors = 0;

exp = [];

for i=0:n_files-1
	load([folder,'\', prefix, num2str(i),'.mat']);
	
	ffModel = FFRandom(2, 2, 2, 10000, result.nt, result.nr); 
	
	total = length(result.timeLine)+1;
	
	P = NPortSourcingProblem(result.timeLine,result.dt,result.chargeData,...
		result.deviceData,result.constraints,ffModel);
	
	ttl = 1000;
	acc = 0;
	while ttl>0
		ttl = ttl-1;
		global lifetime_global;
		lifetime_global=0;
		[~, ~] = P.solve();
		acc = acc + lifetime_global;
		disp(['acc:',num2str(acc),' ttl:',num2str(ttl)]);
	end
	acc
	error('stop');
	exp = [exp, acc/(total*1000)];
end

histogram(exp);