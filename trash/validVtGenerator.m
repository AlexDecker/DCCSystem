%verify if the method of generating valid Vt vectors from img(I) is valid
while true
	for nt = 1:10
		for nr = 1:10
			
			n = nt+nr;
			
			%generating Z components
			w = 2*pi*(100000 + 1000000*rand);%angular frequency
			M = zeros(n);%inductance
			for i=1:n
				for j=i+1:n
					M(i,j) = rand/1000000;
				end
			end
			M = M + M.';%guarantee of symmetry
			R = diag(5*rand(1,n));%diagonal matrix for the resistance values
			Z = R - (1j)*w*M;
				
			Rl = diag([zeros(1,nt),40*rand(1,nr)]);%load resistance
	
			%generating a valir V vector
			V = [50*rand(nt,1);zeros(nr,1)];
	
			%getting the corresponding I vector
			I = (Z+Rl)\V;
	
			%verifying the conditions for the already known valid vector
			%first condition: (equation D)
			error1 = (real(I)-1/w*(M\(R+Rl))*imag(I));
			if(mean(abs(error1))>0.0001)
				disp(' ');disp(' ');
				disp('I found significative error for equation D');
				M
				R
				Rl
				w
				Z
				error1
				error('finishing...');
			end
			
			%the second condition:(equation E)
			error2 = (((1/w*(R(nt+1:end,:)+Rl(nt+1:end,:))/M)*(R+Rl)+w*M(nt+1:end,:))*imag(I));
			if(mean(abs(error2))>0.0001)
				disp(' ');disp(' ');
				disp('I found significative error for equation E');
				M
				R
				Rl
				w
				Z
				error2
				error('finishing...');
			end
			disp(['Error for equation D :',num2str(mean(error1))]);
			disp(['Error for equation E :',num2str(mean(error2))]);
		end
	end
end
