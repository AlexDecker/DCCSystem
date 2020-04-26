clear all;

err=1e-10;

disp('wait...');
for k=1:100000
    [rlTable, convTable, chargeTable] = randomLookupTables();
    dev = DeviceData(rlTable, convTable, chargeTable);
    if ~dev.check()
        error('Invalid device object!!');
    end

    out = rand*convTable(end,2);
    [min_in, max_in] = dev.iConvACDC(out);
	r = rand; in = min_in*r+max_in*(1-r);
    if abs(out - dev.convACDC(in))>err
        error('Error converting ACDC');
    end

    ic = chargeTable(1,2)+rand*(chargeTable(end,2)-chargeTable(1,2));
    [min_input, max_input] = dev.iEffectiveChargeCurrent(ic);
	r = rand; input = min_input*r+max_input*(1-r);
    if abs(ic - dev.effectiveChargeCurrent(input))>err
        error('Error converting the effective charge current');
    end
	
	[min_ir, max_ir] = dev.domain_convACDC();
	r = rand; ir = min_ir*r+max_ir*(1-r);
	in = dev.convACDC(ir);
	ic = dev.effectiveChargeCurrent(in);
	[in0, in1] = dev.iEffectiveChargeCurrent(ic);
	if in+err<in0 || in>in1+err
		error('iEffectiveChargeCurrent');
	end
	[ir0, ~] = dev.iConvACDC(in0);
	[~, ir1] = dev.iConvACDC(in1);
	if ir+err<ir0 || ir>ir1+err
		error('iConvACDC');
	end

	
    %[min_current, max_current] = dev.domain_convACDC();
    %creating a subdomain
    %min_current_sub = min_current + rand*(max_current - min_current);
    %max_current_sub = min_current_sub + rand*(max_current - min_current_sub);
    %creating a valid output variation
    %[min_output, max_output] = dev.domain_iConvACDC();
    %d_dc = min_output + rand^2*(max_output - min_output - 1e-6) + 1e-6;
    %if abs(input1-input2)<di, the difference between the images is limited to d_dc
    %dx = dev.max_variation_convACDC(min_current_sub,max_current_sub,d_dc) - 1e-6;
    %x0 = min_current_sub;
    %y0 = dev.convACDC(x0);
    %while x0+dx <= max_current_sub
    %    y1 = dev.convACDC(x0+dx);
    %    if abs(y0-y1)>d_dc
    %        error('Error with max_variation_convACDC');
    %    end
    %    y0 = y1;
    %    x0 = x0+dx;
    %end

    if mod(k,5000)==0
        fprintf('|');
    end
end
fprintf('\n');
disp('Test performed with SUCCESS!');
    
