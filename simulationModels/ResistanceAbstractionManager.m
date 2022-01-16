% We are testing if each energy consumer connected to a receiving circuit can
% be abstracted by a single resistor whose resistance depends only on the consuming
% current and the charge of the battery.
% You must instantiate a manager and then provide the samples of the problem
% parameters. You must not change the batteries while providing the samples. You must
% change the resistances, inductances, capacitances, voltages, states-of-charge and 
% consuming currents instead.

%CONSIDER A HOMOGENEOUS SYSTEM REGARDING THE BATTERIES!

classdef ResistanceAbstractionManager
	properties
		w % angular frequency of the signals
		max_model_error
		data % Each row: [SOC, ic, ZL]
	end
	methods
		function obj = ResistanceAbstractionManager(f)
			assert(f>0);
			obj.w = 2*pi*f;
			obj.max_model_error = 0;
			obj.data = zeros(0,3);
		end
		
		% Let v be the transmitting voltages (real vector), it be the transmitting
		% currents (phasor, complex vector), ir be the receiving currents (phasor,
		% complex vector), SOC be the vector of States-of-charges (0-1) and ic be
		% the consumed currents (DC, real vector) and the rest of the parameters
		% such that [v;zeros] = [ZT, M; M.', ZR + diag(ZL)][it;ir], where ZL are
		% the equivalent impedances of the consumers.
		function [obj, ZL] = addSample(obj, ZT, ZR, M, v, it, ir, SOC, ic)
			% validating the model
			model_error = abs(v - (ZT * it + M * ir));
			disp('Model error:');
			disp(model_error);
			obj.max_model_error = max(obj.max_model_error, max(model_error));
			disp(['Max model error: ', num2str(obj.max_model_error)]);
			
			ZL = (-M.'*it -ZR*ir)./ir;
			
			assert(length(ic)==length(SOC));
			assert(length(ic)==length(ZL));
			for (i = 1:length(ic))
				obj.data = [obj.data; [ZL(i), ic(i), SOC(i)]];
			end
		end
		
		% Each row of the wave is an observation of the waveform ([time, value]).
		% The function returns a complex number of the phasor representation.
		% It considers the errors in the observations to be iid and normally
		% distributed.
		function [phasor, offset] = inferWaveParameters(obj, wave)
			
			t = wave(:,1);
			y = wave(:,2);
			
			%y = D * b, b = [offset; real(phasor); img(phasor)]
			D = [ones(length(t),1), sin(obj.w * t), cos(obj.w * t)];
			
			%D.'*y = D.'*D * b, so:
			b = (D.'*D)\D.'*y;
			
			offset = b(1);
			phasor = b(2) + 1i*b(3);
		end
	end
	
	methods(Static)
		function test()
			%%%%%% Testing the phasor inference %%%%%%
			
			% parameters
			if (rand < 0.5)
				offset_signal = -1;
			else
				offset_signal = +1;
			end
			offset = offset_signal * (0.5 + rand);
			a = 0.5 + rand;
			f = 1000;
			phi = 2*pi*rand - pi;
			nsamples = 100;
			t = 0;
			for (i = 1:nsamples-1)
				t = [t; rand + t(end)];
			end
			t = 3 * t / (t(end) * f);
			e = normrnd(0, a/5, nsamples, 1);
			y = offset + a*sin(2*pi*f*t + phi) + e;
			
			obj = ResistanceAbstractionManager(f);
			[phasor, offset_] = obj.inferWaveParameters([t, y]);
			a_ = abs(phasor);
			phi_ = atan2(imag(phasor)/a_, real(phasor)/a_);
			y_ = offset_ + a_*sin(2*pi*f*t + phi_);
			
			figure;
			hold on;
			plot(t,y,'o');
			plot(t,y_,'-');
			
			disp(['offset(ref, est): ', num2str(offset), ', ', num2str(offset_)]);
			disp(['amplitude(ref, est): ', num2str(a), ', ', num2str(a_)]);
			disp(['initial phase(ref, est): ', num2str(phi), ', ', num2str(phi_)]);
			
			%%%%%% Testing the ZR inference %%%%%%
			nt = randi(5);
			nr = randi(5);
			ZT = diag(rand(nt,1)) - 1i*rand(nt);
			ZR = diag(rand(nr,1)) - 1i*rand(nr);
			ZL = rand(nr,1) + 1i*rand(nr,1);
			M = -1i*rand(nt,nr);
			Z = [ZT, M; M.', ZR + diag(ZL)];
			v = 10 * rand(nt, 1) - 5;
			i = Z\[v;zeros(nr,1)];
			it = i(1:nt);
			ir = i(nt+1:end);
			[obj, ZL_] = obj.addSample(ZT, ZR, M, v, it, ir, rand(nr,1), rand(nr,1));
			disp('Consumer impedance (reference):');
			disp(ZL);
			disp('Consumer impedance (estimated):');
			disp(ZL_);
		end
	end
end	