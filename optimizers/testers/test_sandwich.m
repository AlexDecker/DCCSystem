clear all;

nSamples = 1;

for i=1:nSamples
	[rlTable,convTable,chargeTable] = randomLookupTables();
	[a1,b1,a2,b2] = sandwich(rlTable(:,1),rlTable(:,2));
	figure;
	hold on;
	plot(rlTable(:,1),rlTable(:,2));
	plot(rlTable(:,1),a1*rlTable(:,1)+b1);
	plot(rlTable(:,1),a2*rlTable(:,1)+b2);
end