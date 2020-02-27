%general arguments
nr = 2;
nt = 2;
nSlots = 2;
dt = 1;
maxV = 15;
sampleSize = 10;%the higher this value, the higher the difficulty

%arguments for the FFDummie
hashSize = 1000;
nSegments = 100;
maxSize = 10000;
thr_top = 10;
thr = 10;
thr_down = 10;
ttl = 100;
ttl_down = 100;

found_solutions = 0;
n_instances = 100;

for i=1:n_instances
    s = rand;%sparsity
    d = rand;%dynamicity
    
    ffModel = FFDummie(hashSize, nSegments, maxSize, thr_top, thr, thr_down, ttl, ttl_down, nt, nr);

    [P, reference_solution] = randomSourcingInstance(nt,nr,nSlots,dt,maxV,sampleSize,s,d,ffModel);
    [result, QLOG_ref] = P.verify(reference_solution);
    if result~=0
        error(['RandomInstance: error number ',num2str(result)]);
    end

    [solveable, solution] = P.solve();

    if solveable
        found_solutions = found_solutions + 1;
        
        [result, Q_LOG] = P.verify(solution);

        if result~=0
            error(['SolveCharging: error number ',num2str(result)]);
        end

        disp(['Instance #',num2str(i),' completed: success'])
    else
        disp(['Instance #',num2str(i),' completed: failure'])
    end
end

disp(['Effectiveness (0-1): ', num2str(found_solutions/n_instances)])
