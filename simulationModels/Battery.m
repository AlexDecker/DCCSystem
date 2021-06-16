classdef Battery
	properties 
		C0
		R0
		Q0
		Q
		C
		R
		j
		i0
		I
		j0
		J
		k
		
		total_charge
		total_power
		open_circuit_voltage
	end
	
	methods
		function obj = Battery()
			obj.C0 = 5.0970e-4;%F
			obj.R0 = 1;%ohm - resistência equivalente do dispositivo
			obj.Q0 = 0;%C

			obj.Q = zeros(11,1);

			obj.C = [2.0466e-3;
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

			obj.R = [1.4024e-3;
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
			% unit vectors (relative positions of the referred variables
			% in the vectorized unified decision variable):
			obj.j = E(1,:);
			obj.i0 = E(2,:);
			obj.I = E(3:12,:);
			obj.j0 = E(13,:);
			obj.J = E(14:23,:);
			obj.k = E(24,:);
			
			% initializing the output variables
			obj.total_charge = 0;
			obj.total_power = 0;
			obj.open_circuit_voltage = 0;
		end
		
		% Updates the battery inner state considering
		% i: average input DC current (A)
		% dt: size of the considered time-slot (s)
		function obj = iterate(obj, i, dt)
			% kirchhoff's current law
			KC = [obj.j + obj.i0;
				 obj.j0 + obj.I(1,:) - obj.i0;
				 obj.J(1,:) + obj.I(2,:) - obj.I(1,:);
				 obj.J(2,:) + obj.I(3,:) - obj.I(2,:);
				 obj.J(3,:) + obj.I(4,:) - obj.I(3,:);
				 obj.J(4,:) + obj.I(5,:) - obj.I(4,:);
				 obj.J(5,:) + obj.I(6,:) - obj.I(5,:);
				 obj.J(6,:) + obj.I(7,:) - obj.I(6,:);
				 obj.J(7,:) + obj.I(8,:) - obj.I(7,:);
				 obj.J(8,:) + obj.I(9,:) - obj.I(8,:);
				 obj.J(9,:) + obj.I(10,:) - obj.I(9,:);
				 obj.J(10,:) + obj.k - obj.I(10,:)];

			KC_y = [i; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

			% kirchoff's voltage law
			KV = [-obj.R0*obj.j + dt/obj.C0*obj.j0;
				  -dt/obj.C0*obj.j0 + obj.R(1)*obj.I(1,:) + dt/obj.C(1)*obj.J(1,:);
				  -dt/obj.C(1)*obj.J(1,:) + obj.R(2)*obj.I(2,:) + dt/obj.C(2)*obj.J(2,:);
				  -dt/obj.C(2)*obj.J(2,:) + obj.R(3)*obj.I(3,:) + dt/obj.C(3)*obj.J(3,:);
				  -dt/obj.C(3)*obj.J(3,:) + obj.R(4)*obj.I(4,:) + dt/obj.C(4)*obj.J(4,:);
				  -dt/obj.C(4)*obj.J(4,:) + obj.R(5)*obj.I(5,:) + dt/obj.C(5)*obj.J(5,:);
				  -dt/obj.C(5)*obj.J(5,:) + obj.R(6)*obj.I(6,:) + dt/obj.C(6)*obj.J(6,:);
				  -dt/obj.C(6)*obj.J(6,:) + obj.R(7)*obj.I(7,:) + dt/obj.C(7)*obj.J(7,:);
				  -dt/obj.C(7)*obj.J(7,:) + obj.R(8)*obj.I(8,:) + dt/obj.C(8)*obj.J(8,:);
				  -dt/obj.C(8)*obj.J(8,:) + obj.R(9)*obj.I(9,:) + dt/obj.C(9)*obj.J(9,:);
				  -dt/obj.C(9)*obj.J(9,:) + obj.R(10)*obj.I(10,:) + dt/obj.C(10)*obj.J(10,:);
				  -dt/obj.C(10)*obj.J(10,:) + obj.R(11)*obj.k + dt/obj.C(11)*obj.k];

			KV_y = [-obj.Q0/obj.C0;
					obj.Q0/obj.C0 - obj.Q(1)/obj.C(1);
					obj.Q(1)/obj.C(1) - obj.Q(2)/obj.C(2);
					obj.Q(2)/obj.C(2) - obj.Q(3)/obj.C(3);
					obj.Q(3)/obj.C(3) - obj.Q(4)/obj.C(4);
					obj.Q(4)/obj.C(4) - obj.Q(5)/obj.C(5);
					obj.Q(5)/obj.C(5) - obj.Q(6)/obj.C(6);
					obj.Q(6)/obj.C(6) - obj.Q(7)/obj.C(7);
					obj.Q(7)/obj.C(7) - obj.Q(8)/obj.C(8);
					obj.Q(8)/obj.C(8) - obj.Q(9)/obj.C(9);
					obj.Q(9)/obj.C(9) - obj.Q(10)/obj.C(10);
					obj.Q(10)/obj.C(10) - obj.Q(11)/obj.C(11)];
					
			K = [KC;KV];
			y = [KC_y;KV_y];

			% values of the currents
			currents = K\y;
			j_ = currents(obj.j==1);
			i0_ = currents(obj.i0==1); 
			I_ = currents(sum(obj.I)==1);
			j0_ = currents(obj.j0==1); 
			J_ = currents(sum(obj.J)==1);
			k_ = currents(obj.k==1);

			% uptading the charges
			obj.Q0 = obj.Q0 + dt*j0_;
			obj.Q = obj.Q + dt*[J_;k_];

			% final parameters
			obj.total_charge = obj.Q0 + sum(obj.Q);
			obj.total_power = obj.R0*j_^2 + sum(obj.R.*[I_;k_].^2) + sum((obj.Q0/obj.C0)*j0_) + sum((obj.Q./obj.C).*[J_;k_]);
			obj.open_circuit_voltage = obj.Q0/obj.C0;
		end
		
		function test(obj, i, dt, T)
			V = zeros(T,1);
			Q = zeros(T,1);
			for t = 1 : T
				disp(obj.open_circuit_voltage)
				disp(V(t))
				V(t) = obj.open_circuit_voltage;
				Q(t) = obj.total_charge;
				obj = obj.iterate(i, dt);
			end
			figure;
			hold on;
			plot((1:T)*dt, V);
			plot((1:T)*dt, Q);
			legend('Voltage (V)', 'Charge (C)');
			xlabel('Time (s)');
			title('Battery evolution');
		end
	end
end