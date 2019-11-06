%This script tests the 'verify' function from NPortChargingProblem.

clear all;

err = 1e-8;%tolerated error

code = 0;%error code (see 'verify' method)

R = diag([1;1;1;5]);

ZT = R(1:2,1:2);
ZR = [0,  -10i;
     -10i, 0]+R(3:4,3:4);
M  = [-10i, -5i;
      -1i, -5i];

Z = [ZT, M.';
     M, ZR];

maxPact = 10;%maximum active power
maxPapp = inf;%maximum apparent power

minCurr = [0;0.875];%min current to charge
maxCurr = [50;50;50;50];%maximum current supported by devices

chargeData.minimum = [0;0];%lower bound
chargeData.initial = [err;err];%in the beginnig
chargeData.threshold = [10;10];%required in the end
chargeData.maximum = [10;10];%complete charge

dt=0.1;

nSlots = 20/dt;

timeLine = [];
for t=1:nSlots
    slot.Z = Z;
    slot.Id = [0;0];
    timeLine = [timeLine;slot];
end

rlTable = [0, 0; 1, chargeData.maximum(1)];
rlCellArray = {rlTable, rlTable};

convCellArray = {[0,0; maxCurr(3),maxCurr(3)],...
                [0,0; minCurr(2)-1e-9,0; minCurr(2),minCurr(2); maxCurr(4),maxCurr(4)]};

NPCP = NPortChargingProblem(timeLine, dt, chargeData, rlCellArray, convCellArray, maxCurr,...
        maxPapp, maxPact);

%GREEDY!!!-----------------------------------------------------------
%%{
q = chargeData.initial;%initial charges
solution = zeros(2,0);

for t=1:nSlots
    disp([num2str(100*t/nSlots),'%']);
    Rl = [interp1(rlCellArray{1}(:,1),rlCellArray{1}(:,2),q(1,end)/chargeData.maximum(1));...
        interp1(rlCellArray{2}(:,1),rlCellArray{2}(:,2),q(2,end)/chargeData.maximum(2))];
    iZ = eye(4)/(timeLine(t).Z+diag([0;0;Rl]));
    c = [-inf;-inf];
    v = zeros(2,1);

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
        %the current vectors adjusted according to the limits
        i0 = i0*min([sqrt(maxPact/(i0'*R*i0)),...
                    sqrt(maxPapp/abs(i0'*v0)),...
                    min(maxCurr./abs(i0))]);
        i1 = i1*min([sqrt(maxPact/(i1'*R*i1)),...
                    sqrt(maxPapp/abs(i1'*v1)),...
                    min(maxCurr./abs(i1))]);
        i2 = i2*min([sqrt(maxPact/(i2'*R*i2)),...
                    sqrt(maxPapp/abs(i2'*v2)),...
                    min(maxCurr./abs(i2))]);
        i3 = i3*min([sqrt(maxPact/(i3'*R*i3)),...
                    sqrt(maxPapp/abs(i3'*v3)),...
                    min(maxCurr./abs(i3))]);
        %the charging currents
        c0 = [interp1(convCellArray{1}(:,1),convCellArray{1}(:,2),abs(i0(3)));...
            interp1(convCellArray{2}(:,1),convCellArray{2}(:,2),abs(i0(4)))]-...
            timeLine(t).Id;
        c1 = [interp1(convCellArray{1}(:,1),convCellArray{1}(:,2),abs(i1(3)));...
            interp1(convCellArray{2}(:,1),convCellArray{2}(:,2),abs(i1(4)))]-...
            timeLine(t).Id;
        c2 = [interp1(convCellArray{1}(:,1),convCellArray{1}(:,2),abs(i2(3)));...
            interp1(convCellArray{2}(:,1),convCellArray{2}(:,2),abs(i2(4)))]-...
            timeLine(t).Id;
        c3 = [interp1(convCellArray{1}(:,1),convCellArray{1}(:,2),abs(i3(3)));...
            interp1(convCellArray{2}(:,1),convCellArray{2}(:,2),abs(i3(4)))]-...
            timeLine(t).Id;
         
        %get the best charging current (only considering the active devices)
        if c0'*(q(:,end)>chargeData.minimum)>sum(c)
            c = c0.*(q(:,end)>chargeData.minimum);
            v = v0;
        end
        if c1'*(q(:,end)>chargeData.minimum)>sum(c)
            c = c1.*(q>chargeData.minimum);
            v = v1;
        end
        if c2'*(q(:,end)>chargeData.minimum)>sum(c)
            c = c2.*(q(:,end)>chargeData.minimum);
            v = v2;
        end
        if c3'*(q(:,end)>chargeData.minimum)>sum(c)
            c = c3.*(q(:,end)>chargeData.minimum);
            v = v3;
        end
    end
    
    %verifying constraints
    if sum(abs(iZ*v)>maxCurr)>0 && code==0
        code = 3;
    end
    if abs((v')*(iZ')*v)>maxPapp && code==0
        code = 4;
    end
    if real((v')*(iZ')*v)>maxPact && code==0
        code = 5;
    end

    solution = [solution, v(1:2)];%log the voltages of this time slot
    q = [q,max(min(q(:,end)+dt*c,chargeData.maximum),chargeData.minimum)];
    
    %charge constraint
    if sum(q(:,end)<=chargeData.minimum)>0 && code==0
        code = 6;
    end
end
%all charges are complete?
if sum(q(:,end)<chargeData.threshold)>0 && code==0
    code = 2;
end

figure;
plot(q');
legend('Device A','Device B');
title('Greedy');
disp(['Code reference: ',num2str(code)]);
disp(['Code answered: ',num2str(verify(NPCP,solution))]);
%}
