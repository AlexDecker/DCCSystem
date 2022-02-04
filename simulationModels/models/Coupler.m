classdef Coupler
	properties
		nt % #transmitters
		nr % #receivers
		mutual_coupler % block
		coupling_helper % mutual inductance generator
	end
	methods
		function coupler = Coupler(systemName, hierarchy, nt, nr, root)
			coupler.nt = nt;
			coupler.nr = nr;
			% local setup (get all coupling references in the pwd)
			coupler.coupling_helper = CouplingHelper([], false, root);
			
			coupler.mutual_coupler.name = [systemName,'/coupler'];

			add_block('powerlib/Elements/Mutual Inductance',...
				coupler.mutual_coupler.name, ...
				'Position', hierarchy.rect ...
			);
			
			set_param(coupler.mutual_coupler.name, 'TypeOfMutual', 'Generalized mutual inductance');
			set_param(coupler.mutual_coupler.name, 'NumberOfWindings', num2str(coupler.nt + coupler.nr));
			
			% for later connection of blocks
			coupler.mutual_coupler.hnd = get_param(coupler.mutual_coupler.name,'PortHandles');
			
			% first set of internal parameters
			coupler = coupler.changeCouplings();
		end
		
		% This function generates and applies a new set of mutual inductances. It also generates
		% the self-inductances and the capacitances in a way to keep the system resonant.
		function coupler = changeCouplings(coupler)
			
			% the mutual induction in a simulated homogeneous system of windings
			M = 4*pi*1e-6 * coupler.coupling_helper.generateMutualInductionMatrix(coupler.nt + coupler.nr);
			
			% the self-induction of each coil
			L = CouplingHelper.referenceSelfInductance();
			
			% the complete inductance matrix
			inductances = -M + L * eye(coupler.nt + coupler.nr);
			
			% setting up the ideal inductors
			set_param(coupler.mutual_coupler.name, 'InductanceMatrix', mat2str(inductances));
			set_param(coupler.mutual_coupler.name, 'ResistanceMatrix', mat2str(zeros(coupler.nt + coupler.nr)));
		end
		
		function hnd = txPositiveHandler(coupler, index)
			hnd = coupler.mutual_coupler.hnd.LConn(index);
		end
		
		function hnd = txNegativeHandler(coupler, index)
			hnd = coupler.mutual_coupler.hnd.RConn(index);
		end
		
		function hnd = rxPositiveHandler(coupler, index)
			hnd = coupler.mutual_coupler.hnd.RConn(coupler.nt + index);
		end
		
		function hnd = rxNegativeHandler(coupler, index)
			hnd = coupler.mutual_coupler.hnd.LConn(coupler.nt + index);
		end
	end
end