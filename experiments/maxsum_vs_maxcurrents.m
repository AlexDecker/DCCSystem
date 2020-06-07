total = 0;
maxsum_success = 0;
maxcurrents_success = 0;

for i=1:350
	load(['baselines\exp_', num2str(i),'.mat']);
	total = total+1;
	maxsum_success = maxsum_success + result.ret(2).success;
	maxcurrents_success = maxcurrents_success + result.ret(1).success;
	
	clear result;
end

p1 = 94/99;%maxsum_success/total
p2 = 45/99;%maxcurrents_success/total

n1 = 99;%total;
n2 = 99;%total;

u = p1-p2;
s = sqrt((p1*(1-p1))/n1 + (p2*(1-p2))/n2);

disp('baselines');
disp(['p-value for maxsum=maxcurrents against maxsum>maxcurrents:', num2str(normcdf(u/s,'upper'))]);