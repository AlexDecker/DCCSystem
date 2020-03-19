%Let Z be the impedance matrix such that [v.',zeros].' = Z [it.', ir.'], where v is the TX voltage vector,
%it is the transmittng current vector and ir is the receiving current vector. This function calculates a
%possible dv such that [(v+dv).',zeros].' = Z [(it+dit).', (ir+dir).'].' and:
% * abs(it+dit)<=It, 
% * abs(ir+dir)=Ir,
% * (v+dv).'*(it+dit)<=P
% * dv real
%This function is based on Newton-Raphson method, which is good since all invelved equations are quadratic
%and the solution dv=0 is near the optimum by hypothesis
function dv_values = fine_adjustment(Z, v, it, ir, It, Ir, P, tolerance, max_iterations)
    
    dv_values = [];
    
    nt = length(it);
    nr = length(ir);

    %mapping the variables into the x vector
    %x = [dv; dit_real; dit_imag; dir_real; dir_imag; power_slack; it_slack]
    
    dv.len = nt;
    dv.begin = 1;
    dv.end = dv.begin + dv.len - 1;
    
    dit_real.len = nt;
    dit_real.begin = dv.end + 1;
    dit_real.end = dit_real.begin + dit_real.len - 1;
    
    dit_imag.len = nt;
    dit_imag.begin = dit_real.end + 1;
    dit_imag.end = dit_imag.begin + dit_imag.len - 1;
    
    dir_real.len = nr;
    dir_real.begin = dit_imag + 1;
    dir_real.end = dir_real.begin + dir_real.len - 1;

    dir_imag.len = nr;
    dir_imag.begin = dir_real.end + 1;
    dir_imag.end = dir_imag.begin + dir_imag.len - 1;

    power_slack.len = 1;
    power_slack.begin = dir_imag.end + 1;
    power_slack.end = power_slack.begin + power_slack.len - 1;

    it_slack.len = nt;
    it_slack.begin = power_slack.end + 1;
    it_slack.end = it_slack.begin + it_slack.len - 1;

    x = zeros(dv.len + dit_real.len + dit_imag.len + dir_real.len + dir_imag.len +...
        power_slack.len + it_slack.len, 1);

    %mapping the equations into F vector
    target_current.len = nr;
    target_current.begin = 1;
    target_current.end = target_current.begin + target_current.len - 1;

    active_power.len = 1;
    active_pover.begin = target_current.end + 1;
    active_power.end = active_power.begin + active_power.len - 1;

    tx_current.len = nt;
    tx_current.begin = active_power.end + 1;
    tx_current.end = tx_current.begin + tx_current.len - 1.

    real_voltage.len = nt+nr;
    real_voltage.begin = tx_current.end + 1;
    real_voltage.end = real_voltage.begin + real_voltage.len - 1;

    imag_voltage.len = nt+nr;
    imag_voltage.begin = real_voltage.end + 1;
    imag_voltage.end = imag_voltage.begin + imag_voltage.len - 1;

    F = zeros(target_current.len + active_power.len + transmitting_current +...
        real_voltage.len + imag_voltage.len, 1);

    %Jacobian
    J = zeros(length(x),length(F));
    
    for k=1:max_iterations
        
        dv_ = x(dc.begin:dv.end);
        
        dit_r = x(dit_real.begin:dit_real.end);
        dit_i = x(dit_imag.begin:dit_imag.end);
        dit = dit_r + 1i*dit_i;

        dir_r = x(dir_real.begin:dir_real.end);
        dir_i = x(dir_imag.begin:dir_imag.end);
        dir = dir_r + 1i*dir_i;

        s_p = x(power_slack.begin:power_slack.end);
        s_t = x(it_slack.begin:it_slack.end);

        %calculating the error vector
        F(target_current.begin:target_current.end) = abs(ir + dir).^2-Ir.^2;
        F(active_power.begin:active_power.end) = (v + dv_).'*(real(it) + dit) + s_p.^2 - P;
        F(tx_current.begin:tx_current.end) = abs(ir + dir).^2 + s_t.^2 - It.^2;
        F(real_voltage.begin:real_voltage.end) = [v + dv_; zeros(nr,1)] - ...
            real(Z)*[real(it)+dit_r ; real(ir)+dir_r] + ...
            imag(Z)*[imag(it)+dit_i ; imag(ir)+dir_i]; 
        F(imag_voltage.begin:imag_voltage.end) = imag(Z) * [real(it)+dit_r ; real(ir)+dir_r] + ...
            real(Z) * [imag(it)+dit_i ; imag(ir)+dir_i];

        %stop condition
        if max(abs(F)) < tolerance
            dv_values = dv_;
            break;
        end

        %calculating the Jacobian
        J(target_current.begin:target_current.end, dv.begin:dv.end) = ...
            zeros(target_current.len,dv.len);
        J(target_current.begin:target_current.end, dit_real.begin:dit_real.end) = ...
            zeros(target_current.len,dit_real.len);
        J(target_current.begin:target_current.end, dit_imag.begin:dit_imag.end) = ...
            zeros(target_current.len,dit_imag.len);
        J(target_current.begin:target_current.end, dir_real.begin:dir_real.end) = ...
            diag(2*real(ir)) + diag(2*dir_r);
        J(target_current.begin:target_current.end, dir_imag.begin:dir_imag.end) = ...
            diag(2*imag(ir)) + diag(2*dir_i);
        J(target_current.begin:target_current.end, power_slack.begin:power_slack.end) = ...
            zeros(target_current.len,power_slack.len);
        J(target_current.begin:target_current.end, it_slack.begin:it_slack.end) = ...
            zeros(target_current.len,it_slack.len);
        
        J(active_power.begin:active_power.end, dv.begin:dv.end) = ...
            diag(dit_r) + diag(real(it));
        J(active_power.begin:active_power.end, dit_real.begin:dit_real.end) = ...
            diag(v) + diag(dv);
        J(active_power.begin:active_power.end, dit_imag.begin:dit_imag.end) = ...
            zeros(active_power.len,dit_imag.len);
        J(active_power.begin:active_power.end, dir_real.begin:dir_real.end) = ...
            zeros(active_power.len,dir_real.len);
        J(active_power.begin:active_power.end, dir_imag.begin:dir_imag.end) = ...
            zeros(active_power.len,dir_imag.len);
        J(active_power.begin:active_power.end, power_slack.begin:power_slack.end) = ...
            diag(2*sp);
        J(active_power.begin:active_power.end, it_slack.begin:it_slack.end) = ...
            zeros(active_power.len,it_slack.len);
        
        J(tx_current.begin:tx_current.end, dv.begin:dv.end) = ...
            zeros(tx_current,dv.len);
        J(tx_current.begin:tx_current.end, dit_real.begin:dit_real.end) = ...
            diag(2*real(it)) + diag(2*dit_r);
        J(tx_current.begin:tx_current.end, dit_imag.begin:dit_imag.end) = ...
            diag(2*imag(it)) + diag(2*dit_i);
        J(tx_current.begin:tx_current.end, dir_real.begin:dir_real.end) = ...
            zeros(tx_current.len,dir_real.len);
        J(tx_current.begin:tx_current.end, dir_imag.begin:dir_imag.end) = ...
            zeros(tx_current.len,dit_imag.len);
        J(tx_current.begin:tx_current.end, power_slack.begin:power_slack.end) = ...
            zeros(tx_current.len,power_slack.len);
        J(tx_current.begin:tx_current.end, it_slack.begin:it_slack.end) = ...
            diag(2*st);

        J(real_voltage.begin:real_voltage.end, dv.begin:dv.end) = ...
            [eye(dv.len);zeros(real_voltage.len,dv.len)]);
        J(real_voltage.begin:real_voltage.end, dit_real.begin:dit_real.end) = ...
            ;
        J(real_voltage.begin:real_voltage.end, dit_imag.begin:dit_imag.end) = ...
            ;
        J(real_voltage.begin:real_voltage.end, dir_real.begin:dir_real.end) = ...
            ;
        J(real_voltage.begin:real_voltage.end, dir_imag.begin:dir_imag.end) = ...
            ;
        J(real_voltage.begin:real_voltage.end, power_slack.begin:power_slack.end) = ...
            zeros(real_voltage.len,power_slack.len);
        J(real_voltage.begin:real_voltage.end, it_slack.begin:it_slack.end) = ...
            zeros(real_voltage.len,it_slack.len);


        %next solution
        x = x - pinv(J)*F;
    end
end
