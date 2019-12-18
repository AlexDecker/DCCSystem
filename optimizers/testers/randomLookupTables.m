%generate random parameters for DeviceData objects
function [rlTable,convTable,chargeTable] = randomLookupTables()
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
        if rand<0.2
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
    while true
        inNeg = chargeTableNeg(1,1)-rand-0.01;
        inPos = chargeTablePos(end,1)+rand+0.01;
        if rand<0.2
            outNeg = chargeTableNeg(1,2);%constant part of the curve
        else
            outNeg = inNeg+rand*(chargeTableNeg(1,2)-inNeg);
        end
        
        chargeTableNeg = [inNeg, outNeg; chargeTableNeg];
        %the input is limited to the maximum converted current
        if inPos>=convTable(end,2) 
            if rand<0.2
                outPos = chargeTablePos(end,2);
            else
                outPos = chargeTablePos(end,2)+...
                    rand*(convTable(end,2)-chargeTablePos(end,2));
            end
            %if the maximum is 0, it may be ommited, since 0 is already guaranteed
            if convTable(end,2)>0
                chargeTablePos = [chargeTablePos; convTable(end,2), outPos];
            end
            break;
        else
            if rand<0.2
                outPos = chargeTablePos(end,2);
            else
                outPos = chargeTablePos(end,2)+rand*(inPos-chargeTablePos(end,2));
            end
            chargeTablePos = [chargeTablePos; inPos, outPos];
        end
        num = num-1;
    end
    chargeTable = [chargeTableNeg; chargeTablePos(2:end,:)];
end
