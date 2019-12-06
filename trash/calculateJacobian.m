function J = calculateJacobian(A, B, C, tolerance, penaltyFactor, w, sb, sc)
    nt = length(B);
    nr = length(A);

    Ja = zeros(nr,3*nt+1);
    for k=1:nr %for each target current constraint
        Ja(k,:) = [w.'*(A{k}+A{k}.'),zeros(1,nt+1)];
    end

    Jb = zeros(nt,3*nt+1);
    for k=1:nt %for each maximum TX current constraint
        Jb(k,:) = [w.'*(B{k}+B{k}.'),zeros(1,nt+1)];
    end

    JC = [w.'*(C+C.'),zeros(1,nt+1)];%power constraint

    %slack constraints
    Js = [zeros(nt,2*nt+1),diag(...
        -penaltyFactor*tolerance*exp(-penaltyFactor*[sb.',sc]))];

    J = [Ja;Jb;Jc;Js];
end
