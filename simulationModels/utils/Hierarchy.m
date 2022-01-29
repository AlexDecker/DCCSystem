% Employed for the hierarchical division of the circuit area
classdef Hierarchy
	properties
		rect
		children
	end
	methods
		function hierarchy = Hierarchy(rect)
			hierarchy.rect = rect;
			hierarchy.children = cell(0);
		end
		
		function hierarchy = horizontalCut(hierarchy, distribution)
			
			if sum(distribution) <= 0 || ~isempty(hierarchy.children)
				error('Invalid cut');
			end
			
			distribution = distribution / sum(distribution);
			
			% top walking dimension
			w = hierarchy.rect(3) - hierarchy.rect(1);
			
			% horizontal cut point
			w_start = hierarchy.rect(1);
			
			for i = 1 : length(distribution)
				
				% horizontal cut point (end)
				w_end = w_start + distribution(i) * w;
				
				hierarchy.children{end + 1} = Hierarchy (...
					[w_start, hierarchy.rect(2), w_end, hierarchy.rect(4)]...
				);
				
				w_start = w_end;
			end
		end
		
		function hierarchy = verticalCut(hierarchy, distribution)
			
			if sum(distribution) <= 0 || ~isempty(hierarchy.children)
				error('Invalid cut');
			end
			
			distribution = distribution / sum(distribution);
			
			% top walking dimension
			h = hierarchy.rect(4) - hierarchy.rect(2);
			
			% vertical cut point
			h_start = hierarchy.rect(2);
			
			for i = 1 : length(distribution)
				
				% vertical cut point (end)
				h_end = h_start + distribution(i) * h;				
				
				% new child-area
				hierarchy.children{end + 1} = Hierarchy(...
					[hierarchy.rect(1), h_start, hierarchy.rect(3), h_end]...
				);
				
				h_start = h_end;
			end
		end
		
		function hierarchy = addPadding(hierarchy, padding)
			
			if ~isempty(hierarchy.children)
				error('Invalid padding insertion');
			end
			
			% top dimensions
			h = hierarchy.rect(3) - hierarchy.rect(1);
			w = hierarchy.rect(4) - hierarchy.rect(2);
			
			hierarchy.children{end + 1} = Hierarchy(...
				hierarchy.rect + [h * padding / 2, w * padding / 2,...
					-h * padding / 2, -w * padding / 2]...
			);
		end
	end
end