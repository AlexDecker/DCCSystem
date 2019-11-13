%This script tests the 'verify' function from NPortChargingProblem.

clear all;

rng('shuffle');

err = 1e-5;%tolerated error
dt=0.25;%integration interval
nSlots = 20/dt;%number of time slots

%generate continuously new test instances
while true
    code = 0;%error code (see 'verify' method)

    maxPact = 10*rand+5;%maximum active power
    maxPapp = maxPact+20*rand;%maximum apparent power

    minCurr = 0.5*rand(2,1);%min current to charge
    maxCurr = 50*rand(4,1)+1;%maximum current supported by devices

    chargeData.minimum = rand(2,1);%lower bound
    chargeData.initial = chargeData.minimum+rand(2,1);%in the beginnig
    chargeData.threshold = chargeData.initial+rand(2,1);%required in the end
    chargeData.maximum = chargeData.threshold+rand(2,1);%complete charge

    timeLine = [];
    for t=1:nSlots
        R = 3*diag(rand(4,1));
        
        MT = rand(2);MR = rand(2);
        
        ZT = R(1:2,1:2)-1i*(MT+MT'-2*diag(diag(MT)))/2;
        ZR = R(3:4,3:4)-1i*(MR+MR'-2*diag(diag(MR)))/2;
        M  = -1i*rand(2);

        Z = [ZT, M.';
             M, ZR];

        slot.Z = Z;
        slot.Id = rand(2,1);
        timeLine = [timeLine;slot];
    end

    %data regarding the load resistance
    rlTable = [0, 0; 1, chargeData.maximum(1)];
    rlCellArray = {rlTable, rlTable};
    
    %current conversion function
    convCellArray = {[0,0; maxCurr(3),maxCurr(3)],...
                    [0,0; minCurr(2)-1e-9,0; minCurr(2),minCurr(2); maxCurr(4),maxCurr(4)]};
    
    NPCP = NPortChargingProblem(timeLine, dt, chargeData, rlCellArray, convCellArray, maxCurr,...
            maxPapp, maxPact, FeasiblePast());
    
    %Starting greedy algorithm
    q = chargeData.initial;%initial charges
    solution = zeros(2,0);
    
    %introducing errors in the bounds
    errPact = (1+(rand<0.01)*rand)*maxPact-100*err;
    errPapp = (1+(rand<0.01)*rand)*maxPapp-100*err;
    errCurr = (1+(rand<0.01)*rand)*maxCurr-100*err;

    for t=1:nSlots
        disp([num2str(100*t/nSlots),'%']);
        Rl = [interp1(rlCellArray{1}(:,1),rlCellArray{1}(:,2),q(1,end)/chargeData.maximum(1));...
            interp1(rlCellArray{2}(:,1),rlCellArray{2}(:,2),q(2,end)/chargeData.maximum(2))];
        iZ = eye(4)/(timeLine(t).Z+diag([0;0;Rl]));
        c = [-inf;-inf];
        v = zeros(4,1);

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
            k0 = min([sqrt(errPact/real(i0'*v0)),...
                        sqrt(errPapp/abs(i0'*v0)),...
                        min(errCurr./abs(i0))]);
            k1 = min([sqrt(errPact/real(i1'*v1)),...
                        sqrt(errPapp/abs(i1'*v1)),...
                        min(errCurr./abs(i1))]);
            k2 = min([sqrt(errPact/real(i2'*v2)),...
                        sqrt(errPapp/abs(i2'*v2)),...
                        min(errCurr./abs(i2))]);
            k3 = min([sqrt(errPact/real(i3'*v3)),...
                        sqrt(errPapp/abs(i3'*v3)),...
                        min(errCurr./abs(i3))]);

            i0 = k0*i0; v0 = k0*v0;
            i1 = k1*i1; v1 = k1*v1;
            i2 = k2*i2; v2 = k2*v3;
            i3 = k3*i3; v3 = k3*v3;

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
                c = c1.*(q(:,end)>chargeData.minimum);
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
        solution = [solution, v(1:2)];%log the voltages of this time slot
        
        %verifying constraints
        if sum(abs(iZ*v)>maxCurr)>0
            code = 3;
            break;
        end
        if abs((v')*(iZ')*v)>maxPapp
            code = 4;
            break;
        end
        if real((v')*(iZ')*v)>maxPact
            code = 5;
            break;
        end
        
        q = [q,max(min(q(:,end)+dt*c,chargeData.maximum),chargeData.minimum)];
        
        %charge constraint
        if sum(q(:,end)<=chargeData.minimum)>0
            code = 6;
            break;
        end
    end
    %all charges are complete?
    if sum(q(:,end)<chargeData.threshold)>0 && code==0
        code = 2;
    end

    [code2, QLog] = verify(NPCP,solution);
    %figure;
    %hold on;
    %plot(QLog.','-');
    %plot(q.','--');
    %legend('result A','result B','reference A','reference B');
    if code~=code2
        error('Incompatible return code');
    else
        disp('Return codes do match!');
    end
end
