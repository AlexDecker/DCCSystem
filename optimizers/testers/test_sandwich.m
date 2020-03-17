clear all;

nSamples = 10000;
maxt = 0;

for i=1:nSamples
	[~,convTable,chargeTable] = randomLookupTables();
    
    tic;
    x = convTable(:,1);
    y = convTable(:,2);
	[a1,b1,a2,b2] = sandwich(x,y);
    if max(a1*x+b1 < y-1e-6 | a2*x+b2 > y+1e-6)~=0
        disp('The bounds are not effective');
        break;
    end
    t = toc;
    maxt = max(maxt, t);

    tic;
    x = chargeTable(:,1);
    y = chargeTable(:,2);
	[a1,b1,a2,b2] = sandwich(x,y);
    if max(a1*x+b1 < y-1e-6 | a2*x+b2 > y+1e-6)~=0
        disp('The bounds are not effective');
        break;
    end
    t = toc;
    maxt = max(maxt, t);
end

disp([num2str(maxt),' seconds at most']);

figure;
hold on;
plot(x,y,'g');
plot(x,a1*x+b1,'b');
plot(x,a2*x+b2,'r');
