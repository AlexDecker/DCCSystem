%this script tests the subSolver, that is, the function for solving the subproblem of
%reachability. Different instances are testes until an error is found.
clear all;

rng('shuffle');

err = 1e-6;%tolerance
penaltyFactor = 5;

maxDim = 10;%maximum number of transmitters or receivers
hardnessLevel = 10;%number of solotions evaluated to get the hardest limits

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

    %testing the returned matrices
    for k=1:1000
        w = 20*(rand(2*nt,1)-0.5);%an ordinary random vector
        x = L*w;%x = [re(it);imag(it)]
        it = x(1:nt)+(1i)*x(nt+1:end);%transmitting current according to the tested function
        ir = H*it;%receiving current according to the tested function
        vt = Z*it;%transmitting voltages according to the tested function

        [e,index] = max(imag(vt));%maximum error for the hypothesis that vt is real
        if e/abs(vt(index))>err %the relative error is more than the tolerated
            error('It returned an out-of-phase voltage vector');
        end

        %now testing if v = ZC*i is respected
        [e,index] = max(abs([vt;zeros(nr,1)]-ZC*[it;ir]));
        if index>nt
            d = 1;
        else
            d = max(vt(index),1);
        end
        if e/d>err
            error('The current values are inconsistent');
        end

        %testing the active power formula
        p = real(it'*vt);
        if p>=1
            if abs(w.'*C*w-p)/p>err
                error('The power value is inconsistent');
            end
        else
            if abs(w.'*C*w-p)>err
                error('The power value is inconsistent');
            end
        end

        %testing the TX current ampĺitude formula
        for l = 1:nt
            if abs(it(l))^2>=1
                if abs(w.'*B{l}*w-abs(it(l))^2)/abs(it(l))^2>err
                    error('TX current inconsistence');
                end
            else
                if abs(w.'*B{l}*w-abs(it(l))^2)>err
                    error('TX current inconsistence');
                end
            end
        end
        
        %testing the RX current ampĺitude formula
        for l = 1:nr
            if abs(ir(l))^2>=1
                if abs(w.'*A{l}*w-abs(ir(l))^2)/abs(ir(l))^2>err
                    error('RX current inconsistence');
                end
            else
                if abs(w.'*A{l}*w-abs(ir(l))^2)>err
                    error('RX current inconsistence');
                end
            end
        end
    end
    
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

    for k=1:1000
        [w,sb,sc] = createInitialSolution(A,B,C,maxIt2,maxP,rand(2*nt,1));
        %verifying the initial solution
        x = L*w;
        it = x(1:nt)+(1i)*x(nt+1:end);
        ir = H*it;
        vt = Z*it;
        
        It2 = abs(it).^2;
        Ir2 = abs(ir).^2;
        p = real(it'*vt);
        
        if maxP>=1
            if (p-maxP)/maxP>err
                error('Initial solution breaks the power constraint');
            end
            if abs(p+sc-maxP)/maxP>err
                error('sc is wrong');
            end
        else
            if p-maxP>err
                error('Initial solution breaks the power constraint');
            end
            if abs(p+sc-maxP)>err
                error('sc is wrong');
            end
        end

        [e,index] = max(It2-maxIt2);
        if maxIt2(index)>=1
            if e/maxIt2(index)>err
                error('Initial solution breaks one of the TX current limits');
            end
            if abs(It2+sb-maxIt2)/maxIt2(index)>err
                error('sb is wrong');
            end
        else
            if e>err
                error('Initial solution breaks one of the TX current limits');
            end
            if abs(It2+sb-maxIt2)>err
                error('sb is wrong');
            end
        end

        f = calculateResidue(A, B, C, err, penaltyFactor, w, sb, sc, absIr2, maxIt2, maxP);

        resIr = Ir2 - absIr2;
        resIt = It2 + sb - maxIt2;
        resP = p + sc - maxP;

        [e,index] = max(abs(f(1:nt+nr+1)-[resIr;resIt;resP]));

        if abs(f(index))>=1
            if e/abs(f(index))>err
                error('Wrong residue');
            end
        else
            if e>err
                error('Wrong residue');
            end
        end

        %testing the slack constraints (<err iff the slack is feasible)
        if sum([sb;sc]>=0 & f(nt+nr+2:end)>err)>0
            error('The slack is ok, but the residue is not ok');
        end
        if sum([sb;sc]<0 & f(nt+nr+2:end)<=err)>0
            error('The slack is not ok, but the residue is ok');
        end

    end
end
