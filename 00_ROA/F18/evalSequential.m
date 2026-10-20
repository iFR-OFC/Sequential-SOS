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
    maxConVio(k) = S.stats.single_iterations{k}.max_con_vio;
    alpha(k) = S.stats.single_iterations{k}.step_length;
    % counterFeasRes(k) = S.stats.single_iterations{k}.counter_restoration_calls;
    % counterSoC(k) = S.stats.single_iterations{k}.counter_correction_calls;
    % counterLineserchSteps(k) = S.stats.single_iterations{k}.counter_linesearch_calls;
end


maxConVio(maxConVio == 0) = eps;

evalSequentialStruct.cost = cost;
evalSequentialStruct.maxConVio = maxConVio;
evalSequentialStruct.alpha = alpha;
evalSequentialStruct.counter_restoration_calls = counterFeasRes;
evalSequentialStruct.counter_correction_calls = counterSoC;
evalSequentialStruct.counter_linesearch_calls = counterLineserchSteps;

end
 