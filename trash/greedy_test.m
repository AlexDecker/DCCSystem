%Counterexample for the greedy charging approach (use always the voltages which maximize
%the sum of the charging currents) that shows sub-optimality.

clear all;

R = diag([1;1;1;3]);

ZT = R(1:2,1:2);
ZR = [0,  -10i;
     -10i, 0]+R(3:4,3:4);
M  = [-10i, -5i;
      -1i, -5i];

Z = [ZT, M.';
     M, ZR];

P = 10;%maximum active power

minI = [0;0.875];%min current to charge
maxQ = [10;10];%complete charge

dt=0.01;

%GREEDY!!!-----------------------------------------------------------
q = zeros(2,1);%initial charges

for t=0:dt:20
    iZ = eye(4)/(Z+diag([0;0;q(:,end)]));
    c = 0;
    i = [];

    %optimizing the sum of the charging currents
    for lambda = linspace(0,1,1000)
        %the norm-1 real voltage vectors
        v0 = [lambda;sqrt(1-lambda^2);0;0];
        v1 = [-lambda;sqrt(1-lambda^2);0;0];
        v2 = [lambda;-sqrt(1-lambda^2);0;0];
        v3 = [-lambda;-sqrt(1-lambda^2);0;0];
        %the corresponding current vectors
        i0 = iZ*v0;
        i1 = iZ*v1;
        i2 = iZ*v2;
        i3 = iZ*v3;
        %the current vectors adjusted according to the active power limit
        i0 = i0*sqrt(P/(i0'*R*i0));
        i1 = i1*sqrt(P/(i1'*R*i1));
        i2 = i2*sqrt(P/(i2'*R*i2));
        i3 = i3*sqrt(P/(i3'*R*i3));
        %the charging currents
        c0 = sum(abs(i0(3:4)));
        c1 = sum(abs(i1(3:4)));
        c2 = sum(abs(i2(3:4)));
        c3 = sum(abs(i3(3:4)));
        if c0>c
            c = c0;
            i = i0;
        end
        if c1>c
            c = c1;
            i = i1;
        end
        if c2>c
            c = c2;
            i = i2;
        end
        if c3>c
            c = c3;
            i = i3;
        end
    end
    q = [q,min(q(:,end)+dt*(abs(i(3:4))>minI).*abs(i(3:4)),maxQ)];
end
figure;
plot(q');
legend('first','last');
title('Greedy');

%DEDICATED----------------------------------------------

q = zeros(2,1);%initial charges

for t=0:dt:5
    iZ = eye(4)/(Z+diag([0;0;q(:,end)]));
    v = [0;1;0;0];
    i = iZ*v;
    i = i*sqrt(P/(i'*R*i));
    q = [q,min(q(:,end)+dt*(abs(i(3:4))>minI).*abs(i(3:4)),maxQ)];
end

for t=6:dt:20
    iZ = eye(4)/(Z+diag([0;0;q(:,end)]));
    v = [1;1;0;0];
    i = iZ*v;
    i = i*sqrt(P/(i'*R*i));
    q = [q,min(q(:,end)+dt*(abs(i(3:4))>minI).*abs(i(3:4)),maxQ)];
end
figure;
plot(q');
legend('first','last');
title('Single Coil');

