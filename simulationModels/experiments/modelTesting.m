% Battery parameters (li-ion)
R_BAT = 0.0018462; %internal battery
K1 = 0.0013733; %polarization constant (V/Ah)
K2 = 0.0013733; %polarization resistance
Q = 5.8783; %maximum battery capacity (Ah)
A = 0.10081; %exponential voltage
B = 9.3941; %exponential capacity (Ah-1)
E0 = 1.3012; %constant voltage

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

%% LEMMA 2: Assuming the voltage drop in LOAD can be written as F + G*ir + H*dir/dt + I*int_{0:t}{ir dtau} + J,
%% with F, G, H, I, J independent of ir, the voltage drop in {D2 + D3} using id as direction reference can be
%% written in the form F' + G'*i + H'*di/dt + I'*int_{0:t}{i dtau} + J', with F', G', H', I', J' independent of {i,ir}

%% LEMMA 3: Filter capacitor, load resistance and battery are also F + G*ir + H*dir/dt + I*int_{0:t}{ir dtau} + J