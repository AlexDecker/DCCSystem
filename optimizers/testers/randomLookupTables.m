%generate random parameters for DeviceData objects
function [rlTable,convTable,chargeTable,maxId,maxIn] = randomLookupTables()
    
    %probability of the next function segment to be constant
    p = 0.2;

    %building a valid RL table
    rlTable = [0, rand];
    while true
        soc = rlTable(end,1)+rand/5+0.01;%do not repeat the previous
        rl = rlTable(end,2)+rand+0.01;%this value also cannot be repeated
        if soc>=1 %the state of charge is limited to 1
            rlTable = [rlTable; 1, rl];
            break;
        end
        rlTable = [rlTable; soc, rl];
    end

    %building a valid conversion table
    convTable = [0, 0];
    num = 10*rand+1;
    while num>0
        in = convTable(end,1)+rand+0.01;
        if rand < p
            out = convTable(end,2);
        else
            out = convTable(end,2)+rand*(in-convTable(end,2));
        end
        convTable = [convTable; in, out];
        num = num-1;
    end

    %building a valid charge table
    chargeTableNeg = [0,0];%the non-positive part
    chargeTablePos = [0,0];%the non-negative part
    outNeg = 0;
    outPos = 0;
    while true
        inNeg = chargeTableNeg(1,1)-rand-0.01;
        inPos = chargeTablePos(end,1)+rand+0.01;
        if rand < p
            outNeg = min(inNeg, outNeg);%constant (if possible) part of the curve
        else
            %otherwise, generates a value that guarantees the curve is monothonically increasing
            %and upper bounded by the input value (you lose more charge than the one provided to
            %the device)
            outNeg = (rand+1)*min(inNeg, outNeg);
        end
        
        chargeTableNeg = [inNeg, outNeg; chargeTableNeg];

        %the input is limited to the maximum converted current
        if inPos >= convTable(end,2) 
            if rand < 1-p
                outPos = outPos + rand*(convTable(end,2) - outPos);
            end
            %if the maximum is 0, it may be ommited, since 0 is already guaranteed
            if convTable(end,2)>0
                chargeTablePos = [chargeTablePos; convTable(end,2), outPos];
            end
            break;
        else
            if rand < 1-p
                outPos = outPos + rand*(inPos - outPos);
            end
            chargeTablePos = [chargeTablePos; inPos, outPos];
        end
        num = num-1;
    end
    chargeTable = [chargeTableNeg; chargeTablePos(2:end,:)];

    %the maximum discharge current supported by the device
    maxId = -chargeTable(1,1);
    %the maximum input current amplitude supported by the device
    maxIn = convTable(end,1);
end
