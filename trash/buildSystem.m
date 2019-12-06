%create the coefficient matrices for the constraints. Using slack variables
%A: target current constraints
%B: maximum transmitting current amplitude constraints
%C: active power constraint
%L: linear transformation matrix for obtaining [Re(it),Im(it)] from w
%Z: linear transformation matrix for obtaining vt from it
%H: linear transformation matrix for obtaining ir from it

function [A, B, C, L, Z, H] = buildSystem(ZT, ZR, M) 
    %some useful definitions
    nt = length(ZT);%number of TX
    nr = length(ZR);%number of RX

    H = -ZR\M;%channel (Ir = H*It)
    Z = ZT+M.'*H;%impedance matrix considering only the TX (v = Z*It)

    %The linear transformation matrix for any 2*nt-real vector w generate a solution which
    %leads to real voltages ([re(It);im(It)]=L*w)
    L = eye(2*nt)-pinv([imag(Z),real(Z)])*[imag(Z),real(Z)];
    
    %Coefficient matrices for evaluations
    %The coefficient matrices for the nr receiving amplitude restrictions
    A = cell(nr);
    for k=1:nr
        hk = H(k,:)'*H(k,:);
        A{k} = L.'*[real(hk), -imag(hk);
                    imag(hk), real(hk)]*L;
    end 

    %The coefficient matrices for the nt transmitting amplitude restrictions
    B = cell(nt);
    for k=1:nt
        bk = [zeros(k-1,1);1;zeros(nt-k,1)];
        B{k} = L.'*diag([bk;bk])*L;
    end 

    %The coefficient matrix for the power constraint
    C = L.'*[real(Z), -imag(Z);
             zeros(nt), zeros(nt)]*L;

end
