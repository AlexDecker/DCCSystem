%This script tests the conjecture that potentially infinit set of feasible futures
%can be represented by a dicrete subset without loss of generality

%THEOREM 1: if the union of the futures that the vectors within the subset can achieve
%is equal to the set of all possible futures for the set, the subset can be used as
%the representant.

clear all;

%creating a single slot considering homogeneous devices------------------------------

nt = 5;nr = 5; %do not forget to test nt>nr, nt=nr and nt<nr
dt = 0.1;%duration of the time slot (s)

%battery capacities
chargeData.maximum = 10*rand(nr,1);

%generating a random generic device
[rlTable, convTable, chargeTable, maxId, maxIn] = randomLookupTables();
dev = DeviceData(rlTable, convTable, chargeTable);

%generating the current consumption of the time slot and the environment parameters
slot = randomTimeLine(nt,nr,1,maxId*ones(nr,1));

%sub-matrices of Z
ZT = slot.Z(1:nt,1:nt);
ZR = slot.Z(nt+1:end,nt+1:end);
M = slot.Z(nt+1:end,1:nt);

%generating a charge vector q0 which achieves q if v is provided-------------------

%Let us suppose a charge vector q0
q0 = rand(nr,1).*chargeData.maximum;

%when the charge vector is q, the rl matrix is
RL = zeros(nr);
for i=1:nr
    RL(i,i) = dev.getRLfromSOC(q0(i)/chargeData.maximum(i));
end

%and a in-phase voltage vector v
v = 30*(rand(nt,1)-0.5);

while true
    %the current vector, in turn, is
    I = (slot.Z+[zeros(nt),zeros(nt,nr);zeros(nr,nt),RL])\[v;zeros(nr,1)];
    it = I(1:nt); ir = I(nt+1:end);
    if abs(ir)<=maxIn
        break;
    else
        v=v/2;
    end
end

%just testing the limiting current

%and, finally, the next charge vector is
q = zeros(nr,1);
for i=1:nr
    input = dev.convACDC(ir(i))-slot.Id(i);
    q(i) = q0(i) + dt*dev.effectiveChargeCurrent(input);
end

%evaluating the formulas for obtaining directly it and ir from v
IT = (ZT-M.'*((ZR+RL)\M))\v;
disp(['Error for equation 1: ', num2str(max(abs(IT-it)))]);
IR = ((M/ZT)*M.'-ZR-RL)\((M/ZT)*v);
disp(['Error for equation 2: ', num2str(max(abs(IR-ir)))]);

%now q0 receives a perturbation dq and the voltage v remains the same-------------
dq = min(0.01*rand(nr,1).*q0,chargeData.maximum-q); %dq is limited to 1% of q

%calculating the variation of the load resistance
dRL = zeros(nr);
for i=1:nr
    dRL(i,i) = dev.getRLfromSOC((q0(i)+dq(i))/chargeData.maximum(i))-RL(i,i);
end
%calculating the variation of the current
I = (slot.Z+[zeros(nt),zeros(nt,nr);zeros(nr,nt),RL+dRL])\[v;zeros(nr,1)];
dit = I(1:nt)-it; dir = I(nt+1:end)-ir;
DIR = ((M/ZT)*M.'-ZR-RL)\dRL*ir;
disp(['Error for equation 3: ', num2str(100*mean(abs(DIR-dir)./abs(dir))),'%']);

%now v also suffers a perturbation dv--------------------------------------------
dv = 0.01*rand(nt,1).*v;
I = (slot.Z+[zeros(nt),zeros(nt,nr);zeros(nr,nt),RL+dRL])\[v+dv;zeros(nr,1)];
dit = I(1:nt)-it; dir = I(nt+1:end)-ir;

B = eye(nr)/((M/ZT)*M.'-ZR-RL);
H = ((M/ZT)*M.'-ZR-RL)\(M/ZT);

DIR = B*dRL*ir+(H+B*dRL*H)*dv;
disp(['Error for equation 4: ', num2str(100*mean(abs(DIR-dir)./abs(dir))),'%']);

%generating the value of dv for q0+dq->q ---------------------------------------


