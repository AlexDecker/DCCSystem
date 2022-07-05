%% Eq to be solved %%
% a(log((ir(t)+i)/(2b)+1) + log((ir(t)-i)/(2b)+1))
%	= A(int{ir(t)}(t) - int{ir(t)}(0)) + B + C*ir(t)

% i = i0 * cos(w*t+phi) -> input current
i0 = 10 * (rand-0.5);
w = 1e+6 * (rand-0.5);
phi = 10 * (rand-0.5);

% load parameters
A = 10 * (rand-0.5);
B = 10 * (rand-0.5);
C = 10 * (rand-0.5);

% diode parameters
a = 10 * (rand-0.5);
b = 10 * (rand-0.5);

% Let ir(t) be the load current.
% Let x(t) = int{ir(t)}(t) - int{ir(t)}(0) = int{ir(t)}(t) - q0
% ir(t) = dx(t)/dt
% dir(t)/dt = dx2(t)/dt

% a(log((ir(t)+i)/(2b)+1) + log((ir(t)-i)/(2b)+1))
%	= A(int{ir(t)}(t) - int{ir(t)}(0)) + B + C*ir(t)

syms x(t)
dx = diff(x);
d2x = diff(dx);
i = i0 * cos(w*t+phi);

ode = a*(log((dx(t)+i)/(2*b)+1) + log((dx(t)-i)/(2*b)+1)) == A*x(t) + B + C*d2x(t);
% int{ir(t)}(0)
%cond1 = x(0) == 0;
% ir(0)
%cond2 = dx(0) == 0;

%conds = [cond1 cond2];
%x_(t) = dsolve(ode,conds);
x_(t) = dsolve(ode);
x_ = simplify(x_)

%https://www.mathworks.com/help/matlab/ref/ode23.html