
%% First part: i(AC current) -> DIODE BRIDGE -> ir -> LOAD
%% _________________
%%      i ->        |
%%                  |
%%                 / \
%%               /    \
%%             D1      D2
%%           /          \
%%   ______/             \_________
%%  | ir ->\             /   ir -> |
%%  |       \          /           |
%%  |        D4      D3            |
%%  |         \    /               |
%%  |          \ /                 |
%%  |           |                  |
%%  |           |                  |
%%__  __________|                  |
%%      <- i                       |
%%  |                              |
%%  |                              |
%%  |____________LOAD______________|
%%
%% The voltage drop in a diode with current x: v = -a*log(x/b+1),
%% where a,b > 0  and x >= 0
%% 
%% LEMMA 1: Assuming that the current in D2 is the same as the current in D4 and
%% the current in D1 is the same as in D3:
%% id = current in D2 and D4
%% id - i = current in D1 and D3
%% ir = id + (id - i) = 2id - i
%% current in D2 and D4 = (ir + i)/2
%% current in D1 and D3 = (ir - i)/2
%% Proof: Let id be the current in D2. The other ones are:
%%        I(D1) = id - i
%%        I(D3) = ir - id
%%        I(D4) = ir - (id - i) = ir + i - id
%%		  From the kirchhoff voltage low (using id as direction reference):
%%		  -a*log(id/b+1) + a*log((ir-id)/b+1) + a*log((ir+i-id)/b+1) - a*log((id-i)/b+1) = 0
%%		  -log((id+b)/b) + log((ir-id+b)/b) + log((ir+i-id+b)/b) - log((id-i+b)/b) = 0 (a is not zero)
%%		  log(((ir-id+b)*(ir+i-id+b)) / ((id+b)*(id-i+b))) = 0
%%		  ((ir-id+b)*(ir+i-id+b)) / ((id+b)*(id-i+b)) = 1
%%		  (ir-id+b)*(ir+i-id+b) = (id+b)*(id-i+b)
%%		  (id-(ir+b))*(id-(ir+i+b)) = (id+b)*(id-(i-b))
%%		  id*(id-(ir+i+b)) - (ir+b)*(id-(ir+i+b)) = id*(id-(i-b)) + b*(id-(i-b))
%%		  id^2 -(ir+i+b)*id - (ir+b)*id + (ir+b)*(ir+i+b) = id^2 - (i-b)*id + b*id - b*(i-b)
%%		  id^2 -i*id -(ir+b)*id - (ir+b)*id + (ir+b)*(ir+i+b) = id^2 - (i-b)*id + b*id - b*(i-b)
%%		  -i*id + (ir+b)*(ir+i+b) = -(i-b)*id + b*id - b*(i-b) (simplifying the terms)
%%		  Spliting id/literals and simplifying the terms:
%%		  (-2*ir-4*b)*id = -2*b*i - ir^2 - ir*i - 2*ir*b
%%		  (-2*ir-4*b)*id = -2*b*(ir+i) - ir*(ir+i)
%%		  -2*(2*b+ir)*id = -(2*b+ir)*(ir+i)
%%		  The term 2*b+ir is guaranteed to be non-zero, since b>0 and ir >= 0. So, we simplify:
%%		  id = (ir + i)/2

%% Testcase:
a = rand + 1;
b = rand + 1;
i = 2 * rand - 1;
ir = rand;
id = (ir + i)/2;
disp('That must be almost zero:');
-a*log(id/b+1) + a*log((ir-id)/b+1) + a*log((ir+i-id)/b+1) - a*log((id-i)/b+1)

clear all

%% LEMMA 2: Assuming the voltage drop in LOAD can be written as F + G*ir + H*dir/dt + I*int_{0:t}{ir dtau} + J,
%% with F, G, H, I, J independent of ir, the voltage drop in {D2 + D3} using id as direction reference can be
%% written in the form F' + G'*i + H'*di/dt + I'*int_{0:t}{i dtau} + J', with F', G', H', I', J' independent of {i,ir}

clear all

%% LEMMA 3: Filter capacitor, load resistance and battery are also F + G*ir + H*dir/dt + I*int_{0:t}{ir dtau} + J
%%
%% ________________________________________BAT INTERNAL RESISTOR_______
%%   ir ->     |     ia ->      |     <- ib                            |
%%             | | ic           | | il                                 |
%%             | V              | V                                    |
%%           FILTER            LOAD                             BATTERY EQUIVALENT
%%         CAPACITOR        RESISTANCE                            VOLTAGE SOURCE
%%             |                |                                      |
%% ____________|________________|______________________________________|
%%
%% Battery parameters (li-ion)
RB = 0.0018462; %internal battery
K1 = 0.0013733; %polarization constant (V/Ah)
K2 = 0.0013733; %polarization resistance
Q = 5.8783; %maximum battery capacity (Ah)
A = 0.10081; %exponential voltage
B = 9.3941; %exponential capacity (Ah-1)
E0 = 1.3012; %constant voltage
%% Charge model: EBAT = E0 - K2*Q / (q + 0.1*Q) * i_filtered - K1*Q / (Q - q) + A*exp(-B*q)
%% Discharge model: EBAT = E0 - K2*Q / (Q - q) * i_filtered - K1*Q / (Q - q) + A*exp(-B*q)
%% For a given q (assuming steady state) and i_filtered = ib
q = rand * Q;
BETA = E0 + A*exp(-B*q) - K1*Q / (Q - q);
ALPHA_C = - K2*Q / (q + 0.1*Q);
ALPHA_D = - K2*Q / (Q - q);
%% Charge model: EBAT = BETA + ALPHA_C * ib
i_filtered = rand;
ib = i_filtered;
disp ('It must be next to zero:');
(BETA + ALPHA_C * ib) - (E0 - K2*Q / (q + 0.1*Q) * i_filtered - K1*Q / (Q - q) + A*exp(-B*q))
%% Discharge model: EBAT = BETA + ALPHA_D * ib
i_filtered = -rand;
ib = i_filtered;
disp ('It must be next to zero:');
(BETA + ALPHA_D * ib) - (E0 - K2*Q / (Q - q) * i_filtered - K1*Q / (Q - q) + A*exp(-B*q))
%%
%% For the set LOAD RESISTANCE + BATTERY
%% The voltage drop considering the direction of ia as reference is equal to:
%% v = RL*il = RB*ib - EBAT
%% From that equality, il = RB/RL*ib - EBAT/RL = RB/RL*ib - (BETA + ALPHA*ib)/RL = (RB - ALPHA)/RL * ib - BETA/RL
%% Using the Kirchhoff current law, il = ia + ib. Thus, ia + ib = (RB - ALPHA)/RL * ib - BETA/RL and
%% ib = -RL/(RL - RB + ALPHA) * ia - BETA/(RL - RB + ALPHA)
%% Then, the voltage drop is
%% v = RL*(ia + ib) = (RL - RL^2/(RL - RB + ALPHA)) * ia - RL * BETA/(RL - RB + ALPHA)
%% Defining F_ = - RL * BETA/(RL - RB + ALPHA) and G_ = (RL - RL^2/(RL - RB + ALPHA))
%% v = F_ + G_*ia
RL = rand;
ia = rand;
ib = -RL/(RL - RB + ALPHA_C) * ia - BETA/(RL - RB + ALPHA_C);
EBAT = BETA + ALPHA_C * ib;
il = ia + ib;
F_ = - RL * BETA/(RL - RB + ALPHA_C);
G_ = RL - RL^2/(RL - RB + ALPHA_C);
disp('Those must be next to zero:');
RL*il - (RB*ib - EBAT)
RL*il - ((RL - RL^2/(RL - RB + ALPHA_C)) * ia - RL * BETA/(RL - RB + ALPHA_C))
RL*il - (F_ + G_*ia)
%%
%% For the set CAPACITOR + LOAD RESISTANCE + BATTERY
%%		Symbol table:
%%			V0: capacitor voltage drop at the beginning of the time slot
%%			C: capacitance
%%			[t0,t1]: time slot
%%			a0, a, b, L: fourier parameters of ia
%%			a0_, a_, b_, L(L is assumed to be the same as ia): fourier parameters of ir
%% The voltage drop considering the direction of ir as reference is equal to:
%% v = -1/C*int{t0:t1}{ic dtau} = F_ + G_*ia
%% Using the Kirchhoff current law, ir = ia + ic <=> ic = ir - ia
%% v = 1/C*int{t0:t1}{ia dtau} - 1/C*int{t0:t1}{ir dtau} = F_ + G_*ia
%% 1/C*int{t0:t1}{ia dtau} - G_*ia = F_ - 1/C*int{t0:t1}{ir dtau}
%% Let the signals be written as fourier series:
%%     ia(tau) = a0/2 + sum{n=1:inf}{a(n)*cos(n*pi/L*tau) + b(n)*sin(n*pi/L*tau)}
%%     ir(tau) = a0_/2 + sum{n=1:inf}{a_(n)*cos(n*pi/L*tau) + b_(n)*sin(n*pi/L*tau)}
%% Whose indefinite integrals are:

%https://www.mathworks.com/help/symbolic/solve-a-single-differential-equation.html

clear all
%%
%% Testing the representation as fourier series:
%% Harmonic form:
w = rand + 0.5 %PRINT
A0 = 2*rand-1 %PRINT
N = 10;
% harmonic coefficients
A = 2*rand(N,1)-1 %PRINT
phi = 2*rand(N,1)-1 %PRINT
% exponential coefficients
Cp = A/2.*exp(-(1i)*phi);
C0 = A0/2;
Cn = flip((Cp.')');
C = [Cn; C0; Cp] %PRINT

M = 1000;
X = linspace(-10*w, 10*w, M).';
S = zeros(M,1);
S_ = zeros(M,1);
for m = 1:M
	x = X(m);
	S(m) = A0/2 + A.'*cos(w*x*(1:N).'-phi);
	S_(m) = C.'*exp(w*(1i)*x*(-N:N).');
end

disp(['Error: ', num2str(max(abs(S-S_)))])

%figure;
%hold on;
%plot(X,S,'-r');
%plot(X,real(S_),':b');
%legend('Harmonic', 'Exponential');

%% Testing the product between two series
Cp2 = rand(N,1)/2.*exp(-(1i)*rand(N,1));
C02 = rand/2;
Cn2 = flip((Cp2.')');
C2 = [Cn2; C02; Cp2] %PRINT
S2 = zeros(M,1);
for m = 1:M
	x = X(m);
	S2(m) = C2.'*exp(w*(1i)*x*(-N:N).');
end

figure;
hold on;
plot(X,real(S_),'--r');
plot(X,real(S2),'--b');
plot(X,real(S_.*S2),'--g');
legend('S', 'S2', 'S*S2');