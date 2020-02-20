%manages load resistance, current conversion and battery charging efficiency
%for a single device. Unlike Device class from the simulator, it does not store
%the state. It is only the encapsulation of lookup tables.
%All following lookup tables must represent monotonically increasing functions.
%Constant intervals are considered to be monotonically increasing (they are also
%monotonically decreasing)
%   *rlTable: Lookup table for load resistance given the SOC
%   *convTable: Lookup table for DC current given the AC input amplitude,
%   from 0A to maxCurr
%   *chargeTable: Lookup table for the effective charge current given the
%   dc input (from the oposite of the maximum discharge current to the maximum
%   converted current)

classdef DeviceData
    properties(Access=public)
		rlTable
		convTable
        chargeTable
    end
    methods
        function obj = DeviceData(rlTable, convTable, chargeTable)
            obj.rlTable = rlTable;
            obj.convTable = convTable;
            obj.chargeTable = chargeTable;
        end
        %returns true if the object is valid
        function b = check(obj)
            b = true;%default
            [n1,n2] = size(obj.rlTable);
            if n1<2 || n2~=2
                disp('Invalid rlTable dimensions');
                b = false;return;
            end
            if obj.rlTable(1,1)~=0 || obj.rlTable(end,1)~=1 ||...
                sum(obj.rlTable(:,2)<0)>0
                disp('Invalid data for rlTable');
                b = false;return;
            end
            for j=2:n1
                if obj.rlTable(j,1)<obj.rlTable(j-1,1) ||...
                    obj.rlTable(j,2)<=obj.rlTable(j-1,2)
                    disp('Invalid data for rlTable');
                    b = false;return;
                end
            end
        
            [n1,n2] = size(obj.convTable);
            if n1<2 || n2~=2
                disp('Invalid convTable dimensions');
                b = false;return;
            end
            %the input aplitude must be between 0 and the real positive maximum.
            %the output must not be greater than the input amplitude.
            if obj.convTable(1,1)~=0 || obj.convTable(end,1)<=0 ||...
                sum(obj.convTable(:,2)>obj.convTable(:,1))>0
                disp('Invalid data for convTable');
                b = false;return;
            end
            for j=2:n1
                if obj.convTable(j,1)<=obj.convTable(j-1,1) ||...
                    obj.convTable(j,2)<obj.convTable(j-1,2)
                    disp('convTable must be monotonically increasing');
                    b = false;return;
                end
            end

            [n1,n2] = size(obj.chargeTable);
            if n1<2 || n2~=2
                disp('Invalid chargeTable dimensions');
                b = false;return;
            end
            %the domain must be from the maximum discharge amplitude to the 
            %maximum charge amplitude. The output must have the same signal
            %as the input (or be zero) and it must be less or equal
            %to the the input (not in absolute value)
            if obj.chargeTable(1,1) >= 0 ||...
                obj.chargeTable(end,1) ~= obj.convTable(end,2) ||...
                sum(obj.chargeTable(:,1) < obj.chargeTable(:,2)) > 0 ||...
                sum(obj.chargeTable(:,2) ~= 0 &...
                    sign(obj.chargeTable(:,1)) ~= sign(obj.chargeTable(:,2))) > 0
                disp('Invalid data for chargeTable');
                b = false;return;
            end
            for j=2:n1
                if obj.chargeTable(j,1) <= obj.chargeTable(j-1,1) ||...
                    obj.chargeTable(j,2) < obj.chargeTable(j-1,2)
                    disp('chargeTable must be monotonically increasing');
                    b = false;return;
                end
            end
        end
        %DOMAIN FUNCTIONS
        %the functions below get lower and upper bounds for the domain of each
        %conversion function
        function [min_SOC, max_SOC] = domain_getRLfromSOC(obj)
            min_SOC = obj.rlTable(1, 1);
            max_SOC = obj.rlTable(end, 1);
        end

        function [min_current, max_current] = domain_convACDC(obj)
            min_current = obj.convTable(1, 1);
            max_current = obj.convTable(end, 1);
        end

        function [min_current, max_current] = domain_effectiveChargeCurrent(obj)
            min_current = obj.chargeTable(1, 1);
            max_current = obj.chargeTable(end, 1);
        end

        function [min_current, max_current] = domain_iConvACDC(obj)
            min_current = obj.convTable(1, 2);
            max_current = obj.convTable(end, 2);
        end

        function [min_current, max_current] = domain_iEffectiveChargeCurrent(obj)
            min_current = obj.chargeTable(1, 2);
            max_current = obj.chargeTable(end, 2);
        end

        %CONVERSION FUNCTIONS
        %get the corresponding load resistance from the actual charge
        function rl = getRLfromSOC(obj,SOC);

            if SOC>1 || SOC<0 || imag(SOC)~=0
                error('getRLfromSOC: SOC is a real number in [0,1]');
            end

            rl = interp1(obj.rlTable(:,1),obj.rlTable(:,2),SOC);
        end

        %get the DC current from an AC input
        function dc_current = convACDC(obj, receiving_current)
            
            %AC->DC without losses
            receiving_current = abs(receiving_current);
            
            %for input validation
            [min_current, max_current] = obj.domain_convACDC;

            if receiving_current > obj.convTable(end,1) || receiving_current < obj.convTable(1,1)
                error(['convACDC: input out of limits. input = ',...
                    num2str(input), '; limits: ', num2str(obj.convTable(1,1)),...
                    ', ', num2str(obj.maxReceivingCurrent())]);
            end

            %losses
            dc_current = interp1(obj.convTable(:,1),obj.convTable(:,2),receiving_current);
        end
        %get the current which will effectively be stored or returned by the
        %battery for a given DC input current (out-discharge_current). Input
        %is required to be more than the discharge current and less than the
        %maximum converted output current
        function ic = effectiveChargeCurrent(obj,input)
            if imag(input)~=0
                error('effectiveChargeCurrent: input must be real');
            end
            if input>obj.chargeTable(end,1) || input<obj.chargeTable(1,1)
                error(['effectiveChargeCurrent: input out of limits. input = ',...
                    num2str(input), '; limits: ', num2str(obj.chargeTable(1,1)),...
                    ', ', num2str(obj.chargeTable(end,1))]);
            end
            ic = interp1(obj.chargeTable(:,1),obj.chargeTable(:,2),input);
        end
        
        %INVERSE CONVERSION FUNCTIONS. Remark: the conversion functions are not
        %always injective, so the inverse functions may not exist.

        %get the input amplitude for a target output (inverse of convACDC)
        function a = iConvACDC(obj, output)
            if output>obj.convTable(end,2) || output<0 || imag(output)~=0
                error(['iConvACDC: invalid parameter output = ',...
                    num2str(output), '; limit: ', num2str(obj.convTable(end,2))]);

            end
            %try to find. Use 'first' if there is more than one occurrence
            i = find(obj.convTable(:,2)==output,1,'first');
            if isempty(i) %did not find, then interpolate
                i0 = find(obj.convTable(:,2)<output,1,'last');
                %i0 cannot be empty, since output is non-negative
                %i1 is i0+1
                a = interp1(obj.convTable(i0:i0+1,2),obj.convTable(i0:i0+1,1), output);
            else
                a = obj.convTable(i,1);%there is an entry with the desired output
            end
        end

        %inverse function of effectiveChargeCurrent
        function input = iEffectiveChargeCurrent(obj, ic)
            if ic>obj.chargeTable(end,2) || ic<obj.chargeTable(1,2) ||...
                imag(ic)~=0
                error(['iEffectiveChargeCurrent: invalid parameter. ic = ',...
                    num2str(ic), '; limits: ', num2str(obj.chargeTable(1,2)),...
                    ', ', num2str(obj.chargeTable(end,2))]);
            end
            %try to find. Use 'first' if there is more than one occurrence
            i = find(obj.chargeTable(:,2)==ic,1,'first');
            if isempty(i) %did not find, then interpolate
                i0 = find(obj.chargeTable(:,2)<ic,1,'last');
                %i0 cannot be empty, since output is non-negative
                %i1 is i0+1
                input = interp1(obj.chargeTable(i0:i0+1,2),obj.chargeTable(i0:i0+1,1), ic);
            else
                input = obj.chargeTable(i,1);%there is an entry with the desired output
            end
        end
    end
end
