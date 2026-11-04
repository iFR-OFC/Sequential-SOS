classdef extractorClass

    properties (SetAccess = private)
        sosProb   % vector of decision-variable polynomials
        sosSol    % solution vector, same order/size as sosProb
    end

    methods
        function obj = extractorClass(sosProb, solx)
    
            obj.sosProb = sosProb(:);
            obj.sosSol= solx(:);
        end

        % get solution 
        function sol = getSol(obj, var)
            
            % get index of decision variable
            idx = obj.findIndex(var);
            sol = reshape(obj.sosSol(idx), size(var));
        end

        function idx = findIndex(obj, var)

            idx = zeros(numel(var), 1);
            for k = 1:numel(var)
                idx(k) = obj.indexOf(var(k));
            end
        end

        function i = indexOf(obj, v)
            % find index of the polynomial in the decision variables
            for i = 1:numel(obj.sosProb)
                if obj.sosProb(i).is_equal(v)
                    return
                end
            end

        end
    end
end
