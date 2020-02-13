%general arguments
nr = 5;
nt = 6;
nSlots = 10;
dt = 1;
maxV = 15;
sampleSize = 100;%the higher this value, the higher the difficulty

%arguments for the FFDummie
hashSize = 1000;
nSegments = 100;
maxSize = 10000;
thr_top = 10;
thr = 10;
thr_down = 10;
ttl = 1000;
ttl_down = 1000;

found_solutions = 0;

for i=1:1
    s = rand;%sparsity
    d = rand;%dynamicity
    
    ffModel = FFDummie(hashSize, nSegments, maxSize, thr_top, thr, thr_down, ttl, ttl_down, nt, nr);

    [P, reference_solution] = randomSourcingInstance(nt,nr,nSlots,dt,maxV,sampleSize,s,d,ffModel);
    [result, ~] = P.verify(reference_solution);
    if result~=0
        error(['RandomInstance: error number ',num2str(result)]);
    end

    [solveable, solution] = P.solve();

    if solveable
        found_solutions = found_solutions + 1;
        
        [result, ~] = P.verify(solution);

        if result~=0
            error(['SolveCharging: error number ',num2str(result)]);
        end

        disp(['Instance #',num2str(i),' completed: success'])
    else
        disp(['Instance #',num2str(i),' completed: failure'])
    end
end

disp(['Effectiveness (0-1): ', num2str(found_solutions/n_instances)])
disp('Finished. SUCCESS');
