%generates a random impedance progress for the system and a random energy consumption.
%nt: number of transmitting devices
%nr: number of receiving devices
%nSlots: number of time slots
%maxId: vector where maxId(i) is the maximum discharge current for device i
%sparsity: a scalar that varies between 0 (very sparce inductance matrix) and +1
%   (very dense)
%dynamicity: a scalar that varies between 0 and +1 and determines if the inductance
%   matrix varies a bit between frames (0) or a lot(+1)

function timeLine = randomTimeLine(nt,nr,nSlots,maxId,sparsity,dynamicity)
    if nt<1 || nr<1 || nSlots<1
        error('randomTimeLine');
    end
    [s1,s2] = size(maxId);
    if sum(imag(maxId)~=0)>0 || s1~=nr || s2~=1 || sum(maxId<0)>0
        error('randomTimeLine: Invalid maxId vector');
    end
    if sparsity<0 || sparsity>1
        error('sparsity must be between -1 and 1');
    end
    if dynamicity<0 || dynamicity>1
        error('dynamicity must be between -1 and 1');
    end

    timeLine = []; 

    %resistance values
    rt = 3*(rand(nt,1)+0.025);
    rr = 5*(rand(nr,1)+0.025);

    %base inductance matrices    
    MT = rand(nt);MR = rand(nr);M = rand(nr,nt);
    %forcing the determined sparsity
    MT = mean(rt)*(1-sparsity)*MT/mean(mean(MT));
    MR = mean(rr)*(1-sparsity)*MR/mean(mean(MR));
    M = mean([rt;rr])*(1-sparsity)*M/mean(mean(M));

    for t=1:nSlots
        %addends for the matrices
        dMT = rand(nt)-0.5;
        dMR = rand(nr)-0.5;
        dM  = rand(nr,nt)-0.5;
        %normalizing the addends
        dMT = (1-dynamicity)*mean(mean(MT))*dMT/mean(mean(dMT));
        dMR = (1-dynamicity)*mean(mean(MR))*dMR/mean(mean(dMR));
        dM  = (1-dynamicity)*mean(mean(M))*dM/mean(mean(dM));
        %base inductance matrices    
        MT = MT+dMT;
        MR = MR+dMR;
        M  = M +dM;
        %forcing the determined sparsity
        MT = mean(rt)*(1-sparsity)*MT/mean(mean(MT));
        MR = mean(rr)*(1-sparsity)*MR/mean(mean(MR));
        M = mean([rt;rr])*(1-sparsity)*M/mean(mean(M));

        %impedance sub-matrices
        %(symmetric, real in the main diagonal and negative imaginary otherwise)
        ZT = diag(rt)-1i*(MT+MT'-2*diag(diag(MT)))/2;
        ZR = diag(rr)-1i*(MR+MR'-2*diag(diag(MR)))/2;
        
        %building the impedance matrix for this time slot
        Z = [ZT,-1i*M.';
             -1i*M, ZR];
        
        slot.Z = Z;
        slot.Id = maxId.*rand(nr,1);%random discharge current, up to 3A
        timeLine = [timeLine;slot];
    end 
end
