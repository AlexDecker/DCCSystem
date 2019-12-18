clear all;

err=1e-6;

disp('wait...');
for k=1:100000
    [rlTable, convTable, chargeTable] = randomLookupTables();
    dev = DeviceData(rlTable, convTable, chargeTable);
    if ~dev.check()
        error('Invalid device object!!');
    end

    out = rand*convTable(end,2);
    in = dev.iConvACDC(out);
    if abs(out- dev.convACDC(in))>err
        error('Error converting ACDC');
    end

    ic = chargeTable(1,2)+rand*(chargeTable(end,2)-chargeTable(1,2));
    input = dev.iEffectiveChargeCurrent(ic);
    if abs(ic - dev.effectiveChargeCurrent(input))>err
        error('Error converting the effective charge current');
    end

    if mod(k,5000)==0
        fprintf('|');
    end
end
fprintf('\n');
disp('Test performed with SUCCESS!');
    
