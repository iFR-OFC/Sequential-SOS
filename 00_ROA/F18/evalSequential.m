function  evalSequentialStruct = evalSequential(S)

% get nmber of entries. Note we also counter the zeroth iteration
nelem = numel(S.stats.single_iterations);

cost      = nan(nelem,1);
maxConVio = nan(nelem,1);
alpha     = nan(nelem,1);
counterSoC = nan(nelem,1);
counterLineserchSteps = nan(nelem,1);
counterFeasRes = nan(nelem,1);


for k = 1:nelem
    cost(k) = S.stats.single_iterations{k}.cost;
    maxConVio(k) = S.stats.single_iterations{k}.maxConVio;
    alpha(k) = S.stats.single_iterations{k}.stepLength;
    % counterFeasRes(k) = S.stats.single_iterations{k}.counterFeasRes;
    % counterSoC(k) = S.stats.single_iterations{k}.counterSoC;
    % counterLineserchSteps(k) = S.stats.single_iterations{k}.counterLineserchSteps;
end


maxConVio(maxConVio == 0) = eps;

evalSequentialStruct.cost = cost;
evalSequentialStruct.maxConVio = maxConVio;
evalSequentialStruct.alpha = alpha;
evalSequentialStruct.counterFeasRes = counterFeasRes;
evalSequentialStruct.counterSoC = counterSoC;
evalSequentialStruct.counterLineserchSteps = counterLineserchSteps;

end
 