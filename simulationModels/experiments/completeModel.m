
%% Circuit scheme:
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
%%  |____________C1________________|
%%  | <- iC1                       | | iH
%%  |                              | V
%%  |                              |
%%  |                              H
%%  |                              |
%%  |                              | 
%%  |____________C2________________|
%%  | <- iC2                       |
%%  |                              | 
%%  |                              | 
%%  |                              | 
%%  |                              | 
%%  |_________-VRL+_______RL_______|
%%  | <- iRL                       |
%%  |                              | 
%%  |                              | 
%%  |                              | 
%%  |                              | 
%%  |_________-VB+________RB_______|
%%    -> iB 
%%
%% Assumed equalities:
%% Shockley diode equation: vD = n * vT * ln {iD/b + 1} where
%%	-> n \in \[1:2\]
%%  -> vT = thermal voltage = k * T / e
%%	   (Boltzman constant, temperature and electron charge)
%%	   for 300 k temperature, vT = 0.0259 V
%%	-> b = bias saturation current, tipically 1e-12 A

% Diode constants
vT = 0.0259;
n = rand + 1;
b = 1e-12;

iD = linspace(1e-10,1,1000); % values <= 0 are undefined or inf
% Dealing with low precision: log(iD + 10^-12) ~ log(iD)
yyaxis left
plot(iD, n * vT * (log(iD) + 12*log(10)) ./ iD);
ylabel('Ohms')
xlabel('A')

yyaxis right
plot(iD, n * vT * (log(iD) + 12*log(10)));
ylabel('Volts')
xlabel('A')

title('Shockley Diode');

figure;

t = linspace(0,10*pi,10000);
iD = abs(cos(t)) + 1e-10; % values <= 0 are undefined or inf
% Dealing with low precision: log(iD + 10^-12) ~ log(iD)
yyaxis left
plot(t, n * vT * (log(iD) + 12*log(10)) ./ iD);
ylabel('Ohms')
xlabel('s')

yyaxis right
plot(t, n * vT * (log(iD) + 12*log(10)));
ylabel('Volts')
xlabel('s')

title('Shockley Diode: AC over time in bridge');


