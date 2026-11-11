%% -----------------------------------------------------------------------
%
%  Short Description: This script generates the evaluation for a desired
%                     benchmarks.  It reads in the .mat file from the
%                     desired benchmark problem. The user can decide which
%                     plots are generated
%
%
%% -----------------------------------------------------------------------


function [] = evaluate_benchmark(matFileName)

load(matFileName)


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
    solverNames{s} = ['Sequential (' seqMethodNames{s} ')'];
end
solverNames{end} = 'Coordinate-descent';

% evaluate 
plotPerformanceProfiles(seqMethodStructs, ...
                        cdSolverStructs, ...
                        solverNames, ...
                        commonOpts.maxSolveTime, ...
                        probNames)


% final results table accumulated over all runs and all problems
resultsTable = cell2table(tableRows,'VariableNames',tableVariableNames);
disp(resultsTable)

%% Plot convergence for each run
for kk = 1:numel(seqSolverStructs_regularization)

if enable_plot_convergence
    if enable_second_sequential
        
        evalSequentialStruct_regularization = evalSequential(seqSolverStructs_regularization{kk});
        evalSequentialStruct_regularization.hess_est_method = hess_est_method_1;

        evalSequentialStruct_hessApprox2                 = evalSequential(seqSolverStructs_hessApprox2{kk});
        evalSequentialStruct_hessApprox2.hess_est_method = hess_est_method_2;



        % plot convergence
        plotConvergence(probNames{kk}, ...
            enable_second_sequential,...
            evalSequentialStruct_regularization,...
            evalSequentialStruct_hessApprox2,...
            cdSolverStructs{kk})
    else

        evalSequentialStruct_regularization = evalSequential(seqSolverStructs_regularization{kk});
        evalSequentialStruct_regularization.hess_est_method = hess_est_method_1;

        % plot convergence
        plotConvergence(probNames{kk}, ...
            enable_second_sequential,...
            evalSequentialStruct_VDP_regularization,...
            [],...
            cdSolverStructs{kk})
    end
end

end

end