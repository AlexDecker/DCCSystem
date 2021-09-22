classdef CouplingHelper
	properties
	
		% empirical cumulative distribution
		cdf_x
		cdf_y
		
		n_coils % number of devices modeled
	end
	methods
		% CouplingHelper generates new inductive coupling matrices based on
		% provided reference data.
		%  references. Cell Array containing the reference matrices
		%  plot_input. If true, plots the correlation matrix of the reference
		%		 matrices and the histogram of the mutual inductances.
		function obj = CouplingHelper(references, plot_input)
			
			if isempty(references)
				references = obj.loadFromDisk();
			end
			
			sample_size = length(references);
			obj.n_coils = length(references{1});
			
			% number of independent (supposed) and identically distributed (due to symmetry)
			% elements per matrix (symmetric and with empty main diagonal)
			iid_elements = (obj.n_coils^2 - obj.n_coils) / 2;
			
			% every independent element in a single vector
			sample_index = 1;
			reference_sample = zeros(sample_size * iid_elements, 1);
			
			% vectorized independent elements, but grouped by origin matrix (in a consistent order).
			% Used for independence testing.
			vec_row = 1;
			vectorized_matrices = zeros(sample_size, iid_elements);
			
			for matrix_index = 1:sample_size
				
				n_coils_ = length(references{matrix_index});
				
				if obj.n_coils ~= n_coils_
					error('Inconsistent references');
				end
				
				% vectorized version of the matrix
				row_index = 1;
				row = zeros(iid_elements, 1);
				
				for i = 1 : obj.n_coils
					for j = i+1 : obj.n_coils
						reference_sample(sample_index) = references{matrix_index}(i,j);
						sample_index = sample_index + 1;
						% vectorized version of the matrix
						row(row_index) = references{matrix_index}(i,j);
						row_index = row_index + 1;
					end
				end
							
				% for sake of independence testing
				vectorized_matrices(vec_row, :) = row.';
				vec_row = vec_row + 1;
			end
			
			% cumulative distribution
			obj.cdf_x = sort(reference_sample);
			obj.cdf_y = linspace(0,1,length(obj.cdf_x)).';
			
			if plot_input
				P = corr(vectorized_matrices,'type','spearman');
				
				figure;
				image(P,'CDataMapping','scaled');
				title('Reference correlationn');
				colorbar;
				
				figure;
				histogram(reference_sample, 'normalization', 'pdf');
				title('PDF');
				
				figure;
				hold on;
				histogram(reference_sample, 'normalization', 'cdf');
				plot(obj.cdf_x, obj.cdf_y);
				title('CDF');
			end
		end
		
		function references = loadFromDisk(obj)
		
			% loading all availiable neumann-based matrices in this directory
			n_instances = 0;
			
			mat = dir('*.mat');
			
			references = cell(length(mat),1);
			
			for q = 1:length(mat)
				disp(['Loading ',mat(q).name, '...']);
				load(mat(q).name);
				references{q} = M;
			end

			if isempty(references)
				error('No .mat file in this directory');
			end
			
		end
		
		% Sample from the target distribution using the inverse-cumulative method
		function m = generateMutualInduction(obj, n_samples)
			% Approximating the inverse via linear interpolation
			m = interp1(obj.cdf_y, obj.cdf_x, rand(n_samples, 1));
		end
		
		% Verify if the generated matrices are right according the references
		function test(obj)
			
			sample_size = length(obj.cdf_x); %cdf_x has the whole original sample
		
			test_sample = obj.generateMutualInduction(sample_size);
			
			% kolmogorov-smirnov test
			[~, p] = kstest2(obj.cdf_x, test_sample);
			
			disp(['p-value (H0: two samples are identically distributed): ', num2str(p)]);
			
			nbins = ceil(sample_size / 10);
			
			figure;
			hold on;
			histogram(obj.cdf_x, nbins, 'normalization', 'pdf');
			histogram(test_sample, nbins, 'normalization', 'pdf');
			legend('reference','test');
		end
		
		function M = generateMutualInductionMatrix(obj, n_coils)
			
			M = zeros(n_coils);
			
			iid_elements = (n_coils^2 - n_coils) / 2;
			
			sample = obj.generateMutualInduction(iid_elements);
			
			k = 1;
			for i = 1 : n_coils
				for j = i+1 : n_coils
					M(i,j) = sample(k);
					M(j,i) = sample(k);
					k = k + 1;
				end
			end
			
		end
		
	end
	
	properties(Constant)
		% coil dimensions (we admit only homogeneous systems)
		R2 = 0.15; % external radius
		R1 = 0.05;% internal radius
		N = 25; % number of turns
		pts = 750; % coil resolution
		mi = 4*pi*10e-7; % vacuum magnectic permeability
		wire_radius = 0.002;
	end
	
	methods(Static)
		% Based on 
		%@article{mohan1999simple,
		%  title={Simple accurate expressions for planar spiral inductances},
		%  author={Mohan, Sunderarajan S and del Mar Hershenson, Maria and Boyd, Stephen P and Lee, Thomas H},
		%  journal={IEEE Journal of solid-state circuits},
		%  volume={34},
		%  number={10},
		%  pages={1419--1424},
		%  year={1999},
		%  publisher={IEEE}
		%}
		% Generates the self inductance of the reference coil (Planar winding) in Henrys
		function L = referenceSelfInductance()
			% average ratio
			davg = (CouplingHelper.R2 + CouplingHelper.R1);
			% fill ratio
			p = (CouplingHelper.R2 - CouplingHelper.R1) / (CouplingHelper.R2 + CouplingHelper.R1);
			% Coefficients for Current Sheet Expression (circular coil)
			c1 = 1.00;
			c2 = 2.46;
			c3 = 0.00;
			c4 = 0.20;
			% less than 8% error for s < 3w:
			L = 0.5 * CouplingHelper.mi * CouplingHelper.N^2 * davg * c1 * (log(c2/p) + c3*p + c4*p^2);
		end
	
		% Creating inductance matrices using the Neumann integral (more reliable)
		function neumannSampling()
			rng('shuffle')

			% sample size (increase this value if you want more instances, use 0 if you already have enough ones)
			sample_size = 0;%100;

			% number of coils
			n_coils = 6;

			prototype = SpiralPlanarCoil(...
				CouplingHelper.R2,CouplingHelper.R1,...
				CouplingHelper.N,CouplingHelper.wire_radius,...
				CouplingHelper.pts...
			);

			% simulation area restricted to a 1 m origin-centered cube
			maxV = 0.5;
			maxR = pi/2;
			
			n = 1;
			while true
				group_list = [];
		
				% creating the coils as planar-spirals
				for c = 1:n_coils
					coil = translateCoil(prototype,unifrnd(-maxV,maxV),...
											unifrnd(-maxV,maxV),unifrnd(-maxV,maxV));
					group.coils.obj = rotateCoilX(rotateCoilY(coil,unifrnd(-maxR,maxR)),unifrnd(-maxR,maxR));
					group.R = -1;group.C = -1;
					group_list = [group_list; group];
					
				end
				
				% creating an environment (encapsulates the impedance matrix and calculates
				% the couplings. Using vacuum magnetic permeability and a dummie frequency
				environment = Environment(group_list, 1e+6, CouplingHelper.mi);
				
				if check(environment)
					% calculate the mutual inductance matriz with no previous information
					environment = environment.evalM(-ones(n_coils) + diag(ones(n_coils,1)));
					
					M = environment.M;
					
					filename = [num2str(round(rand*1000000000)),'.mat'];
					save(filename,'M');
					
					disp(['Instance number ',num2str(n),' saved successfully as ', filename]);
				end
				
				n = n + 1;
			end
		end
	end
end