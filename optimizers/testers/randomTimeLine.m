%generates a random impedance progress for the system and a random energy consumption.
function timeLine = randomTimeLine(nt,nr,nSlots,maxId)
    if nt<1 || nr<1 || nSlots<1
        error('randomTimeLine');
    end
    [s1,s2] = size(maxId);
    if sum(imag(maxId)~=0)>0 || s1~=nr || s2~=1 || sum(maxId<0)>0
        error('randomTimeLine: Invalid maxId vector');
    end
    timeLine = []; 
    for t=1:nSlots
        %resistance values
        RT = 3*diag(rand(nt,1)+0.025);
        RR = 5*diag(rand(nr,1)+0.025);

        %base inductance matrices    
        MT = rand(nt);MR = rand(nr);

        %impedance sub-matrices
        %(symmetric, real in the main diagonal and negative imaginary otherwise)
        ZT = RT-1i*(MT+MT'-2*diag(diag(MT)))/2;
        ZR = RR-1i*(MR+MR'-2*diag(diag(MR)))/2;
        M  = -1i*rand(nr,nt);
        
        %building the impedance matrix for this time slot
        Z = [ZT, M.';
             M, ZR];
        
        slot.Z = Z;
        slot.Id = maxId.*rand(nr,1);%random discharge current, up to 3A
        timeLine = [timeLine;slot];
    end 
end
