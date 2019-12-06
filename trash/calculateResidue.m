%if the values are all less than the threshold, the solution satisfies all constraints

function f = calculateResidue(A, B, C, tolerance, penaltyFactor, w, sb, sc, absIr2, maxIt2, maxP)
    nt = length(B);
    nr = length(A);

    fa = zeros(nr,1);
    for k=1:nr %target current constraints
        fa(k) = w.'*A{k}*w - absIr2(k);
    end
    
    fb = zeros(nt,1);
    for k=1:nt %maximum transmitting current constraints
        fb(k) = w.'*B{k}*w + sb(k) - maxIt2(k);
    end

    fc = w.'*C*w + sc - maxP;%power constraint

    %slack constraints
    fsb = tolerance*exp(-penaltyFactor*sb);
    fsc = tolerance*exp(-penaltyFactor*sc);

    f = [fa;fb;fc;fsb;fsc];
end
