clear all;
rng('shuffle');

err = 1e-9;

n = 6;%half transmitters, half receivers

aux = rand(n/2);
aux = aux-diag(diag(aux));
ZT = diag(10*rand(1,n/2))-(1i)*5*(aux+aux');

aux = rand(n/2);
aux = aux-diag(diag(aux));
ZR = diag(10*rand(1,n/2))-(1i)*5*(aux+aux');

Mrt = -(1i)*rand(n/2);

%the restricted problem
if true
    ZT = real(ZT) + (1i)*(imag(Mrt).'/imag(ZR))*imag(Mrt);
    ZR = (1i)*imag(ZR);
end

Z = [ZT, Mrt.';
    Mrt, ZR];

Vt = 10*(rand(n/2,1)-0.5);

I = Z\[Vt;zeros(n/2,1)];
It = I(1:n/2);
Ir = I(n/2+1:end);

if sum(abs(Vt-(Mrt.'-(ZT/Mrt)*ZR)*Ir))/sum(abs(Vt))>err
    error('equation 1');
end

T = real(ZT);Q = imag(ZT);
R = real(ZR);S = imag(ZR);
m = eye(n/2)/imag(Mrt);

if sum(abs(Vt-((1i)*((eye(n/2)/m).'+T*m*R-Q*m*S)-Q*m*R-T*m*S)*Ir))/sum(abs(Vt))>err
    error('equation 2');
end

H = -Mrt\ZR;
if sum(abs(It-H*Ir))>err
    error('equation 3');
end

H = -m*S;
%if the conditional of the restricted problem is false, this one will probably yield an error
if sum(abs(It-H*Ir))>err
    error('equation 4');
end
