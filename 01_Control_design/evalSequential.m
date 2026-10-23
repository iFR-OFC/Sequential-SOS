function  evalSequentialStruct = evalSequential(S)

% get number of entries. Note we also counter the zeroth iteration
nelem = numel(S.stats.single_iterations);

cost                      = nan(nelem,1);
maxConVio                 = nan(nelem,1);
alpha                     = nan(nelem,1);
counter_correction_calls  = nan(nelem,1);
counter_linesearch_calls  = nan(nelem,1);
counter_restoration_calls = nan(nelem,1);


for k = 1:nelem
    cost(k)                         = S.stats.single_iterations{k}.cost;
    maxConVio(k)                    = S.stats.single_iterations{k}.max_con_vio;
    alpha(k)                        = S.stats.single_iterations{k}.step_length;
    counter_restoration_calls(k)    = S.stats.single_iterations{k}.counter_restoration_calls;
    counter_correction_calls(k)     = S.stats.single_iterations{k}.counter_correction_calls;
    counter_linesearch_calls(k)     = S.stats.single_iterations{k}.counter_linesearch_calls;
end


maxConVio(maxConVio == 0) = eps;

evalSequentialStruct.cost                      = cost;
evalSequentialStruct.maxConVio                 = maxConVio;
evalSequentialStruct.alpha                     = alpha;
evalSequentialStruct.counter_restoration_calls = counter_restoration_calls;
evalSequentialStruct.counter_correction_calls  = counter_correction_calls;
evalSequentialStruct.counter_linesearch_calls  = counter_linesearch_calls;
evalSequentialStruct.iterations                = S.stats.iterations;
evalSequentialStruct.options                   = S.stats.options;
end
 