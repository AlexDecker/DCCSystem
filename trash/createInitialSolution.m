%creates a good initial solution w for which the maximum amplitude and maximum power constraints
%are respected. sb is a vector of slack variables for the B constraints and sc is a slack variable
%for C constraint. w0 is a base vector, in such a way that w will be calculated as a multiple of w0
function [w,sb,sc] = createInitialSolution(A,B,C,maxIt2,maxP,w0)
    nt = length(B);
    nr = length(A);
    %what is the maximum multiplier m1 for w=m1*w0 to respect each amplitude constraint?
    m1 = zeros(nt,1);
    for k=1:nt
        %(m1*w0).'*B*(m1*w0)=maxIt^2 <=> m1^2 = maxIt^2/w0'Bw0
        m1(k) = sqrt(maxIt2(k)/(w0.'*B{k}*w0));
    end
    %what is the maximum multiplier m2 for w=m2*w0 to respect the power constraint?
    %(m2*w0).'*C*(m2*w0)=maxP <=> m1^2 = maxP/w0'Cw0
    m2 = sqrt(maxP/(w0.'*C*w0));
    %getting a multiplier which guarantees both constraint sets to be satistied
    m = min(min(m1),m2);
    %getting a good initial solution
    w = rand*m*w0;
    %calculating the slack variables
    sb = zeros(nt,1);
    for k=1:nt
        sb(k) = maxIt2(k) - w.'*B{k}*w;
    end
    sc = maxP - w.'*C*w;
end
