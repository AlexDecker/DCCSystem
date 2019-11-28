%This function solves a sub-problem of both N-Port Charging Problem and N-Port Power Sourcing
%Problem. Given a reference for the amplitude of the receiving currents, return a transmitting
%real voltage vector which respects the active power and the transmitting current constraints
%and leads to the desired receiving current.
%The function returns an empty vector if no valid solution was found.
%For the naming conventions, see the reference material
%ttl1: maximum number of attempts considering different initial solutions
%ttl2: maximum number of iterations per attempt
function v = subSolver(ZT, ZR, M, maxIt, maxP, absIr, ttl1, ttl2)

    %some useful definitions
    nt = length(ZT);%number of TX
    nr = length(ZR);%number of RX
    H = -ZR\M;%channel (Ir = H*It)
    Z = ZT+M.'*H;%impedance matrix considering only the TX (v = Z*It)
    %The linear transformation matrix for any 2*nt-real vector w generate a solution which
    %leads to real voltages ([re(It);im(It)]=L*w)
    L = eye(2*nt)-pinv([imag(Z),real(Z)])*[imag(Z),real(Z)];
    
    %Coefficient matrices for evaluations
    %The coefficient matrices for the nr receiving amplitude restrictions
    A = cell(nr);
    for k=1:nr
        hk = H(k,:)'*H(k,:);
        A{k} = L.'*[real(hk), -imag(hk);
                    imag(hk), real(hk)]*L;
    end
    %The coefficient matrices for the nt transmitting amplitude restrictions
    B = cell(nt);
    for k=1:nt
        bk = [zeros(k-1,1);1;zeros(nt-k,1)];
        B{k} = diag([bk;bk]);
    end
    %The coefficient matrix for the power constraint
    C = L.'*[real(Z), -imag(Z);
             zeros(nt), zeros(nt)]*L;

    while true
        %create an initial solution. It must respect the maximum TX amplitude and power
        %constraints.
        w0 = rand(2*nt,1); %any solution
        %what is the maximum multiplier m1 for w=m1*w0 to respect each amplitude constraint?
        m1 = zeros(nt,1);
        for k=1:nt
            %(m1*w0).'*B*(m1*w0)=maxIt^2 <=> m1^2 = maxIt^2/w0'Bw0
            m1(k) = sqrt(maxIt(k)^2/(w0.'*B{k}*w0));
        end
        %what is the maximum multiplier m2 for w=m2*w0 to respect the power constraint?
        %(m2*w0).'*C*(m2*w0)=maxP <=> m1^2 = maxP/w0'Cw0
        m2 = sqrt(maxP/(w0.'*C*w0));
        %getting a multiplier which guarantees both constraint sets to be satistied
        m = min(min(m1),m2);
        %getting a good initial solution
        w = m*w0;
        %calculating the slack variables
        fb = zeros(nt,1);
        for k=1:nt
            fb(k) = maxIt(k)^2 - w.'*B{k}*w;
        end
        fc = maxP - w.'*C*w;
        
        %restarting time to leave
        ttl = ttl2;

        %error threshold for failure so far
        merr = inf;

        while true
            %calculating the residues (nt+nr+1 for the regular constraints and
            %nt+1 for the slack ones
            f = zeros(2*nt+nr+2);
            for i=1:nr
                f(k) = w.'*A{k}*w - absIr(k)^2;
            end
            for i=1:nt
                f(k+nr) = w.'*B{k}*w + fb(k) - maxIt(k)^2;
            end
            f(nt+nr+1) = w.'*C*w + fc - maxP;
            %VARIAVEIS DE FOLGA AGR

            if ttl==1 || err>=merr
                w = [];
                break;
            end
            ttl = ttl-1;
        end

        if ttl1==1 %no more attempts
            break; %finish with all you have
        end
        ttl1 = ttl1-1; %less one attempt...
    end
end
