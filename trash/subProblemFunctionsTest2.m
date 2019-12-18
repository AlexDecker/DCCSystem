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
	
	%random voltage vector
    vt = 30*(rand(nt,1)-0.5);

    %getting the accurate currents
    I = ZC\[vt;zeros(nr,1)];
    It = I(1:nt);
    Ir = I(nt+1:end);

    %the decision variable
    x = [real(It);imag(It)];

    H = -ZR\M;%channel (Ir = H*It)    
    for k=1:nr
        hk = H(k,:);
        Ir2 = x.'*[real(hk'*hk), -imag(hk'*hk);
                    imag(hk'*hk), real(hk'*hk)]*x;
        if abs(Ir(k))^2>=1
            if abs(Ir2-abs(Ir(k))^2)/abs(Ir(k))^2>err
                error('First test failed');
            end
        else
            if abs(Ir2-abs(Ir(k))^2)>err
                error('First test failed');
            end
        end

        if sum(sum(abs(real(hk).'*real(hk)+imag(hk).'*imag(hk) - real(hk'*hk))>err))>0
            error('Second test failed');
        end

        if sum(sum(abs(real(hk).'*imag(hk)-imag(hk).'*real(hk) - imag(hk'*hk))>err))>0
            error('Third test failed');
        end
        
        Ir2 = x.'*[real(hk).';-imag(hk).']*[real(hk),-imag(hk)]*x + ...
            x.'*[imag(hk).';real(hk).']*[imag(hk),real(hk)]*x;
        if abs(Ir(k))^2>=1
            if abs(Ir2-abs(Ir(k))^2)/abs(Ir(k))^2>err
                error('Fourth test failed');
            end
        else
            if abs(Ir2-abs(Ir(k))^2)>err
                error('Fourth test failed');
            end
        end
        
        k1 = -1/sqrt(2)*abs(Ir(k));
        k2 = 1/sqrt(2)*abs(Ir(k));
        l1 = -1/sqrt(2)*abs(Ir(k));
        l2 = 1/sqrt(2)*abs(Ir(k));

        eq = [x;1].'*[real(hk).';-imag(hk).';k1]*[real(hk),-imag(hk),k2]*[x;1] + ...
            [x;1].'*[imag(hk).';real(hk).';l1]*[imag(hk),real(hk),l2]*[x;1];

        if abs(eq)>err
            error('Literal absorption test failed');
        end


    end
end
