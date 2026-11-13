%% -----------------------------------------------------------------------
%
%  Short Description: This script generates the evaluation for multiple
%                     desired benchmarks.  It reads in the .mat files
%                     of the desired benchmark problems.
%
%
%% -----------------------------------------------------------------------


function [] = evaluate_benchmark_multiple(matFileNames)

numBenchmarks = numel(matFileNames);

seqMethodStructsAll     = cell(1, numBenchmarks);
cdSolverStructsAll      = cell(1, numBenchmarks);
solverNamesAll          = cell(1, numBenchmarks);
benchmarkNamesAll       = cell(1, numBenchmarks);
solversAll              = cell(1, numBenchmarks);
maxSolveTimeAll         = cell(1, numBenchmarks);
probNamesAll            = cell(1, numBenchmarks);
tableRowsAll            = cell(1, numBenchmarks);
tableVariableNamesAll   = cell(1, numBenchmarks);

for mm = 1:numBenchmarks

    load(matFileNames{mm})

    fprintf('\n')
    fprintf('-----------------------------------------------------------------------------------------------------\n')
    fprintf('Benchmark %d of %d: %s (solver: %s)\n', mm, numBenchmarks, name_of_benchmark, commonOpts.solver)
    fprintf('Date of benchmark: %s\n', date_of_benchmark)
    fprintf('MATLAB version: %s\n', matlabVersion)
    fprintf('CasADi version: %s\n', casadiVersion)
    fprintf('-----------------------------------------------------------------------------------------------------\n')
    fprintf('\n')

% sequential runs stored per Hessian estimation method 
seqMethodStructs = {seqSolverStructs_regularization};
seqMethodNames   = {hess_est_method_1};
if enable_second_sequential
    seqMethodStructs{end+1} = seqSolverStructs_hessApprox2;
    seqMethodNames{end+1}   = hess_est_method_2;
end

% solver names: one sequential solver per Hessian estimation method plus coordinate descent
solverNames = cell(1, numel(seqMethodNames) + 1);
for s = 1:numel(seqMethodNames)
    solverNames{s} = ['Sequential (' seqMethodNames{s} ', ' commonOpts.solver ')'];
end
solverNames{end} = ['Coordinate-descent',' (' commonOpts.solver ')'];

seqMethodStructsAll{mm}     = seqMethodStructs;
cdSolverStructsAll{mm}      = cdSolverStructs;
solverNamesAll{mm}          = solverNames;
benchmarkNamesAll{mm}       = name_of_benchmark;
solversAll{mm}              = commonOpts.solver;
maxSolveTimeAll{mm}         = commonOpts.maxSolveTime;
probNamesAll{mm}            = probNames;
tableRowsAll{mm}            = tableRows;
tableVariableNamesAll{mm}   = tableVariableNames;

%% Plot convergence for each run
for kk = 1:numel(seqSolverStructs_regularization)

if enable_plot_convergence
    if enable_second_sequential
        
        evalSequentialStruct_regularization = evalSequential(seqSolverStructs_regularization{kk});
        evalSequentialStruct_regularization.hess_est_method = hess_est_method_1;

        evalSequentialStruct_hessApprox2                 = evalSequential(seqSolverStructs_hessApprox2{kk});
        evalSequentialStruct_hessApprox2.hess_est_method = hess_est_method_2;



        % plot convergence
        % plotConvergence(probNames{kk}, ...
        %     enable_second_sequential,...
        %     evalSequentialStruct_regularization,...
        %     evalSequentialStruct_hessApprox2,...
        %     cdSolverStructs{kk})
    else

        evalSequentialStruct_regularization = evalSequential(seqSolverStructs_regularization{kk});
        evalSequentialStruct_regularization.hess_est_method = hess_est_method_1;

        % plot convergence
        % plotConvergence(probNames{kk}, ...
        %     enable_second_sequential,...
        %     evalSequentialStruct_VDP_regularization,...
        %     [],...
        %     cdSolverStructs{kk})
    end
end

end

end

% evaluate 
plotPerformanceProfiles_multiple(seqMethodStructsAll, ...
                                 cdSolverStructsAll, ...
                                 solverNamesAll, ...
                                 benchmarkNamesAll, ...
                                 maxSolveTimeAll, ...
                                 probNamesAll)


% final results table accumulated over all runs and all problems
for mm = 1:numBenchmarks

    fprintf('\n')
    fprintf('-----------------------------------------------------------------------------------------------------\n')
    fprintf('%s (solver: %s)\n', benchmarkNamesAll{mm}, solversAll{mm})
    fprintf('-----------------------------------------------------------------------------------------------------\n')
    resultsTable = cell2table(tableRowsAll{mm},'VariableNames',tableVariableNamesAll{mm});
    disp(resultsTable)

end

end
