clear all;

nt = randi(10);
nr = randi(10);
timeLine = randomTimeLine(nt,nr,1,zeros(nr,1),rand,rand);
Z = timeLine.Z;
ZT = Z(1:nt,1:nt);
M = Z(nt+1:end,1:nt);
ZR = Z(nt+1:end,nt+1:end);
V = 15*rand(nt,1);
I = Z\[V;zeros(nr,1)];
IT = I(1:nt);
IR = I(nt+1:end);

%The expansion of the effect of noise for current amplitude
dIR = 0.01*rand(nr,1).*IR; % noise of 1% magnitude
re_dIR = (abs(IR+dIR).^2 - abs(IR).^2 - 2*imag(dIR).*imag(IR))./(2*real(IR))
real(dIR)

