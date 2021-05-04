C0 = 5.0970e-4;%F
R0 = 1;%ohm - resistência equivalente do dispositivo
Q0 = 0;%C
dt = 1;%s
i = 1; %A, input current

Q = zeros(11,1);

C = [2.0466e-3;
	 3.3532e-3;
	 5.8817e-3;
	 1.0305e-2;
	 1.7529e-2;
	 2.8886e-2;
	 4.6367e-2;
	 7.3346e-2;
	 1.1701e-1;
	 1.9984e-1;
	 4.9492e-1];

R = [1.4024e-3;
	 2.6077e-3;
	 4.4226e-3;
	 7.8087e-3;
	 1.3495e-2;
	 2.2587e-2;
	 3.6697e-2;
	 5.8362e-2;
	 9.2332e-2;
	 1.5055e-1;
	 2.8428e-1];

% variables: [j, I[0..10], J[0..10], k]
E = eye(24);
% unit vectors:
j = E(1,:);
i0 = E(2,:);
I = E(3:12,:);
j0 = E(13,:);
J = E(14:23,:);
k = E(24,:);

% kirchhoff's current law
KC = [j + i0;
	 j0 + I(1,:) - i0;
	 J(1,:) + I(2,:) - I(1,:);
	 J(2,:) + I(3,:) - I(2,:);
	 J(3,:) + I(4,:) - I(3,:);
	 J(4,:) + I(5,:) - I(4,:);
	 J(5,:) + I(6,:) - I(5,:);
	 J(6,:) + I(7,:) - I(6,:);
	 J(7,:) + I(8,:) - I(7,:);
	 J(8,:) + I(9,:) - I(8,:);
	 J(9,:) + I(10,:) - I(9,:);
	 J(10,:) + k - I(10,:)];

KC_y = [i; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

% kirchoff's voltage law
KV = [-R0*j + dt/C0*j0;
	  -dt/C0*j0 + R(1)*I(1,:) + dt/C(1)*J(1,:);
	  -dt/C(1)*J(1,:) + R(2)*I(2,:) + dt/C(2)*J(2,:);
	  -dt/C(2)*J(2,:) + R(3)*I(3,:) + dt/C(3)*J(3,:);
	  -dt/C(3)*J(3,:) + R(4)*I(4,:) + dt/C(4)*J(4,:);
	  -dt/C(4)*J(4,:) + R(5)*I(5,:) + dt/C(5)*J(5,:);
	  -dt/C(5)*J(5,:) + R(6)*I(6,:) + dt/C(6)*J(6,:);
	  -dt/C(6)*J(6,:) + R(7)*I(7,:) + dt/C(7)*J(7,:);
	  -dt/C(7)*J(7,:) + R(8)*I(8,:) + dt/C(8)*J(8,:);
	  -dt/C(8)*J(8,:) + R(9)*I(9,:) + dt/C(9)*J(9,:);
	  -dt/C(9)*J(9,:) + R(10)*I(10,:) + dt/C(10)*J(10,:);
	  -dt/C(10)*J(10,:) + R(11)*k + dt/C(11)*k];

KV_y = [-Q0/C0;
		Q0/C0 - Q(1)/C(1);
		Q(1)/C(1) - Q(2)/C(2);
		Q(2)/C(2) - Q(3)/C(3);
		Q(3)/C(3) - Q(4)/C(4);
		Q(4)/C(4) - Q(5)/C(5);
		Q(5)/C(5) - Q(6)/C(6);
		Q(6)/C(6) - Q(7)/C(7);
		Q(7)/C(7) - Q(8)/C(8);
		Q(8)/C(8) - Q(9)/C(9);
		Q(9)/C(9) - Q(10)/C(10);
		Q(10)/C(10) - Q(11)/C(11)];
		
K = [KC;KV];
y = [KC_y;KV_y];

% values of the currents
currents = K\y;
j_ = currents(j==1);
i0_ = currents(i0==1); 
I_ = currents(sum(I)==1);
j0_ = currents(j0==1); 
J_ = currents(sum(J)==1);
k_ = currents(k==1);

% uptading the charges
Q0 = Q0 + dt*j0_;
Q = Q + dt*[J_;k_];

% final parameters
total_charge = Q0 + sum(Q);
total_power = R0*j_^2 + sum(R.*[I_;k_].^2) + sum((Q0/C0)*j0) + sum((Q./C).*[J_;k_]);
open_circuit_voltage = Q0/C0;