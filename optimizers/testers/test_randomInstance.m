nr = 5;
nt = 6;
nSlots = 10;
dt = 1;
maxV = 15;
sampleSize = 100;%the higher this value, the higher the difficulty

for i=1:1000
    s = rand;%sparsity
    d = rand;%dynamicity
    
    ffModel = FeasibleFuture();    
    [P, sol] = randomSourcingInstance(nt,nr,nSlots,dt,maxV,sampleSize,s,d,ffModel);
    [result, ~] = P.verify(sol);

    if mod(i,10)==0
        disp([num2str(i/10),'%']);
    end

    if result~=0
        error(['Error number ',num2str(result)]);
    end
end
disp('Finished. SUCCESS');
