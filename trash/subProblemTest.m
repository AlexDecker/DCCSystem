%this script tests the subSolver, that is, the function for solving the subproblem of
%reachability. Different instances are testes until an error is found.
clear all;

rng('shuffle');

err = 1e-6;%tolerance
penaltyFactor = 5;

maxDim = 10;%maximum number of transmitters or receivers
hardnessLevel = 10;%number of solotions evaluated to get the hardest limits

%building the log matrix (grouped using the number of TX and RX as index)
for nr = 1:maxDim
    for nt = 1:maxDim
        LOG(nr,nt).timeList = [];%execution times
        LOG(nr,nt).iterationList = [];%number of iterations
        LOG(nr,nt).converged = 0;%number of intances which converged
        LOG(nr,nt).diverged = 0;%number of instances which diverged from the result
    end
end

num = 1;

while true
    %generating valid numbers of TX and RX
    nr = round((maxDim-1)*rand)+1;
    nt = round((maxDim-1)*rand)+1;

    disp(['Starting execution number ',num2str(num),...
        ' with ',num2str(nt),' TX and ',num2str(nr),' RX']);
    num = num+1;

    RT = 3*diag(rand(nt,1))+0.01;%TX resistances (between 0.01 ohm and 3.01 ohms)
    RR = 23*diag(rand(nr,1))+0.01;%RX resistances (between 0.01 ohm and 3.01 ohms with
    %load resistance between 0 and 20 ohms -- see magmimo paper)
    
    MTb = rand(nt);MRb = rand(nr);%matrices used as base to build valid  TX-TX and RX-RX couplings
    ZT = RT-1i*(MTb+MTb'-2*diag(diag(MTb)))/2;%TX impedance
    ZR = RR-1i*(MRb+MRb'-2*diag(diag(MRb)))/2;%RX impedance
    M  = -1i*rand(nr,nt);%TX-RX coupling
    
    %complete impedance matrix (v = ZC*i)
    ZC = [ZT, M.';
          M, ZR];
    
    %the first subSolver sub-function to test: coefficient generator
    [A, B, C, L, Z, H] =  buildSystem(ZT,ZR,M);
    
    maxIt2 = inf*ones(nt,1);
    absIr2 = zeros(nr,1);
    maxP = inf;
    hardest = 0;
    %creating the limits
    for k=1:hardnessLevel
        w = 20*(rand(2*nt,1)-0.5);%an ordinary random vector
        x = L*w;
        it = x(1:nt)+(1i)*x(nt+1:end);
        ir = H*it;
        vt = Z*it;
        
        It2 = abs(it).^2;
        Ir2 = abs(ir).^2;
        p = real(it'*vt);

        %hardness score (something like the efficiency)
        h = mean(Ir2)/mean(It2);
        if h>hardest
            maxIt2 = It2;
            maxIr2 = Ir2;
            maxP = p;
        end
    end
    
    itt1 = 100;
    success = false;
    while itt1>0
        itt2 = 100;
        [w,sb,sc] = createInitialSolution(A,B,C,maxIt2,maxP,rand(2*nt,1)-0.5);
        while itt2>0
            f = calculateResidue(A, B, C, err, penaltyFactor, w, sb, sc, absIr2, maxIt2, maxP);
            
            %evaluating the error
            [e, index] = max(abs(f));
            if f(index)>=1
                if e/abs(f(index))<err
                    success = true;
                    break;
                end
            else
                if e<err
                    success = true;
                    break;
                end
            end
            
            %if the solution is not ready, iterate.
            J = calculateJacobian(A, B, C, err, penaltyFactor, w, sb, sc);

            %verify if the Jacobian has inf or NaN
            if sum(sum(isnan(J) | J==inf | J==-inf))>0
                break;%failure
            end

            x = [w;sb;sc] - pinv(J)*f;
            w = x(1:2*nt);
            sb = x(2*nt+1:end-1);
            sc = x(end);

            itt2 = itt2-1;
        end
        if success
            break;
        end
        itt1 = itt1-1;
    end

    LOG(nr,nt).timeList = [LOG(nr,nt).timeList,toc];
    LOG(nr,nt).iterationList = [LOG(nr,nt).iterationList,100*(100-itt1)+100-itt2];
    if success
        LOG(nr,nt).converged = LOG(nr,nt).converged+1;
        disp('SUCCESS');
    else
        LOG(nr,nt).diverged = LOG(nr,nt).diverged+1;
        disp('FAILURE');
    end
end
