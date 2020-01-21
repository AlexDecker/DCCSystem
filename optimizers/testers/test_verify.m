%This script tests the 'verify' function from NPortChargingProblem.
clear all;

err = 1e-6;%tolerated error
dt=0.25;%integration interval
nSlots = 20/dt;%number of time slots

for k=1:100
    code = 0;%error code (see 'verify' method)

    maxPact = 10*rand+5;%maximum active power [5,15]
    maxPapp = maxPact+20*rand;%maximum apparent power [5,35]


    chargeData.minimum = rand(2,1);%lower bound [0,1]
    chargeData.initial = chargeData.minimum+rand(2,1);%in the beginnig [0,2]
    chargeData.threshold = chargeData.initial+rand(2,1);%required in the end [0,3]
    chargeData.maximum = chargeData.threshold+rand(2,1);%complete charge [0,4]

    [rlTable1, convTable1, chargeTable1, maxId1, maxIn1] = randomLookupTables();                            
    [rlTable2, convTable2, chargeTable2, maxId2, maxIn2] = randomLookupTables();                            
    dev1 = DeviceData(rlTable1, convTable1, chargeTable1);
    dev2 = DeviceData(rlTable2, convTable2, chargeTable2);
    
    s = rand;%sparsity
    d = rand;%dynamicity
    timeLine = randomTimeLine(2,2,nSlots,[maxId1;maxId2],s,d);

    minCurr = 0.5*rand(2,1);%min current to charge [0,0.5]
    maxCurr = [50*rand(2,1)+1;maxIn1;maxIn2];%maximum current supported by devices [1,50]

    constraints.maxCurr = maxCurr;
    constraints.maxPapp = maxPapp;
    constraints.maxPact = maxPact;

    NPCP = NPortPowerProblems(timeLine, dt, chargeData, [dev1;dev2], constraints, FeasibleFuture());
    
    %Starting greedy algorithm
    q = chargeData.initial;%initial charges
    solution = zeros(2,0);
    
    %introducing errors in the bounds
    errPact = (1+(rand<0.01)*rand)*maxPact-100*err;
    errPapp = (1+(rand<0.01)*rand)*maxPapp-100*err;
    errCurr = (1+(rand<0.01)*rand)*maxCurr-100*err;

    for t=1:nSlots
        disp([num2str(100*t/nSlots),'%']);
        Rl = [interp1(rlTable1(:,1),rlTable1(:,2),q(1,end)/chargeData.maximum(1));...
            interp1(rlTable2(:,1),rlTable2(:,2),q(2,end)/chargeData.maximum(2))];
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
            c0 = [interp1(convTable1(:,1),convTable1(:,2),abs(i0(3)));...
                interp1(convTable2(:,1),convTable2(:,2),abs(i0(4)))]-...
                timeLine(t).Id;
            c1 = [interp1(convTable1(:,1),convTable1(:,2),abs(i1(3)));...
                interp1(convTable2(:,1),convTable2(:,2),abs(i1(4)))]-...
                timeLine(t).Id;
            c2 = [interp1(convTable1(:,1),convTable1(:,2),abs(i2(3)));...
                interp1(convTable2(:,1),convTable2(:,2),abs(i2(4)))]-...
                timeLine(t).Id;
            c3 = [interp1(convTable1(:,1),convTable1(:,2),abs(i3(3)));...
                interp1(convTable2(:,1),convTable2(:,2),abs(i3(4)))]-...
                timeLine(t).Id;

            %the effective charging currents
            c0 = [interp1(chargeTable1(:,1),chargeTable1(:,2),c0(1));
                interp1(chargeTable2(:,1),chargeTable2(:,2),c0(2))];
            c1 = [interp1(chargeTable1(:,1),chargeTable1(:,2),c1(1));
                interp1(chargeTable2(:,1),chargeTable2(:,2),c1(2))];
            c2 = [interp1(chargeTable1(:,1),chargeTable1(:,2),c2(1));
                interp1(chargeTable2(:,1),chargeTable2(:,2),c2(2))];
            c3 = [interp1(chargeTable1(:,1),chargeTable1(:,2),c3(1));
                interp1(chargeTable2(:,1),chargeTable2(:,2),c3(2))];

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
    %{
    figure;
    hold on;
    plot(QLog.','-');
    plot(q.','--');
    plot((chargeData.minimum.*ones(2,length(q))).','.r');
    plot((chargeData.maximum.*ones(2,length(q))).','.g');
    legend('result A','result B','reference A','reference B');
    %}
    disp(['reference code: ',num2str(code),', result code: ',num2str(code2)]);
    if code~=code2
        error('Incompatible return code');
    end
end
