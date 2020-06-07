total = 0;
pareto_success = 0;
maxsum_success = 0;

for i=0:99
	load(['power_sourcing\exp_', num2str(i),'.mat']);
	
	%if result.nSegments > 13
		total = total+1;
		pareto_success = pareto_success + result.ret(3).success;
		maxsum_success = maxsum_success + result.ret(4).success;
	%end
	
	clear result;
end

p1 = pareto_success/total
p2 = maxsum_success/total

n1 = total;
n2 = total;

u = p1-p2;
s = sqrt((p1*(1-p1))/n1 + (p2*(1-p2))/n2);

disp('power_sourcing');
disp(['p-value for pareto=maxsum against pareto>maxsum:', num2str(normcdf(u/s,'upper'))]);

total = 0;
pareto_success = 0;
maxsum_success = 0;

for i=1:127
	load(['final\exp_', num2str(i),'.mat']);
	
	if result.nSegments > 15
		total = total+1;
		pareto_success = pareto_success + result.ret(1).success;
		maxsum_success = maxsum_success + result.ret(2).success;
	end
	
	clear result;
end

p1 = pareto_success/total
p2 = maxsum_success/total

n1 = total;
n2 = total;

u = p1-p2;
s = sqrt((p1*(1-p1))/n1 + (p2*(1-p2))/n2);

disp('final');
disp(['p-value for pareto=maxsum against pareto>maxsum:', num2str(normcdf(u/s,'upper'))]);