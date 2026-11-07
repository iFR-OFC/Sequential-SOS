%% ***********************************************************************
%
%  Barrier Certificate Benchmark Test Runner
%
%  Benchmark 1: Van-der-Pol Oscillator
%
%  Note 1: Each runs uses various combinations of polynomial degrees. For SOS
%  problem to be feasible or to improve solution quality higher degree are
%  needed. Hence, with a higher degree more decision variables occur and
%  the Gram matrices become larger. Therefore, we can check the behavior
%  for over several benchmark problems, i.e., we have sub-benchmarks so to
%  say that increase our benchmark set.
%
%  Note 2: The user might select a second hessian approximation method for
%  sequential. This is optional, but interesting in terms of convergence
%  speed.
%
%
% ***********************************************************************
set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');

close all
clear
clc


maxNumCompThreads(1);

% add the problem specific folders and the shared helper functions
runnerDir = pwd;

addPaths  = { ...
    fullfile(runnerDir,'VDP'), ...
    fullfile(runnerDir,'..','helperFunctions')};

for ii = 1:numel(addPaths)
    addpath(addPaths{ii})
end

% enable plotting of convergence for each run
enable_plot_convergence  = false;
enable_second_sequential = false;

% Hessian estimation methods of the sequential solver to compare
hess_est_method_1 = 'regularization'; % our default (see also sequentialCommon.m)
hess_est_method_2 = 'mirroring';      % 'BFGS', 'levenberg', 'mirroring'

% common options shared by the methods
commonOpts.maxiter       = 350;
commonOpts.tolerance_opt = 1e-4;
commonOpts.maxSolveTime  = 3600;
commonOpts.solver        = 'mosek';
commonOpts.miniter       = 5; % minimum iterations before an aborted CD run still counts as successful

% MOSEK settings
commonOpts.sdp_solver_opts.mosek_param = struct( 'MSK_IPAR_NUM_THREADS', 1); 
commonOpts.sdp_solver_opts.mosek_echo  = 0;


% variable names of the results table (one row per run)
tableVariableNames = {'Problem Name','no. con. (SDP)','no. dec. (SDP)','Size Gram matrices', ...
                        'Build Time (CD)','Build Time (SQ)','Solve Time (CD)','Solve Time (SQ)', ...
                        'Iterations (CD)','Iterations (SQ)','Final Cost (CD)','Final Cost (SQ)', ...
                        'Solved (CD)','Solved (SQ)'};

% extend the results table by the sequential run with the second Hessian
% estimation method (if enable_second_sequential == 1)
if enable_second_sequential
    tableVariableNames = [tableVariableNames, ...
        ['Build Time (SQ ' hess_est_method_2 ')'], ...
        ['Solve Time (SQ ' hess_est_method_2 ')'], ...
        ['Iterations (SQ ' hess_est_method_2 ')'], ...
        ['Final Cost (SQ ' hess_est_method_2 ')'], ...
        ['Solved (SQ ' hess_est_method_2 ')']];
end

% initialize the cell array that accumulates one row per run of the
% results table
tableRows = cell(0, numel(tableVariableNames));

% initialize the cell arrays of stored runs (one entry per run)
seqSolverStructs_regularization = {};
seqSolverStructs_hessApprox2    = {};
cdSolverStructs                 = {};
probNames                       = {}; % problem name per run (for the evaluation output)


%% VDP

% define the problem instances based on different polynomial degrees
degV = {0:2 , 0:4};
degK = {0:1 , 0:4};

for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_VDP = generateVDP_probStruct(degV{kk},degK{kk});
    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Reachbility Analysis Benchmark for %s Instance %d of %d. \n', probStruct_VDP.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_VDP(probStruct_VDP);

    % solve VDP example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_VDP_regularization = solve_reachability_seq_bench_fun(probStruct_VDP,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_VDP_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_VDP_regularization                 = evalSequential(solStruct_seq_VDP_regularization);
    evalSequentialStruct_VDP_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_VDP_regularization.probNumber      = kk;

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_VDP_hessApprox2 = solve_reachability_seq_bench_fun(probStruct_VDP,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_VDP_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_VDP_hessApprox2 = evalSequential(solStruct_seq_VDP_hessApprox2);
        evalSequentialStruct_VDP_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_VDP_hessApprox2.probNumber = kk;
    end

    % solve VDP example with two states using custom coordinate-descent
    solStruct_CD_VDP = solve_reachability_CD_bench_fun(probStruct_VDP,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_VDP;
    probNames{end+1}       = sprintf('%s (instance %d of %d)', probStruct_VDP.Name, kk, numel(degV));

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows,probStruct_VDP,solStruct_CD_VDP, ...
            solStruct_seq_VDP_regularization,solStruct_seq_VDP_hessApprox2,true);
    else
        tableRows = appendResultRow(tableRows,probStruct_VDP,solStruct_CD_VDP, ...
                                    solStruct_seq_VDP_regularization,[],true);
    end

    if enable_plot_convergence
        if enable_second_sequential
        % Evaluate
        plotConvergence(probStruct_VDP, ...
                        enable_second_sequential,...
                        evalSequentialStruct_VDP_regularization,...
                        evalSequentialStruct_VDP_hessApprox2,...
                        solStruct_CD_VDP)
        else
            % Evaluate
            plotConvergence(probStruct_VDP, ...
                enable_second_sequential,...
                evalSequentialStruct_VDP_regularization,...
                [],...
                solStruct_CD_VDP)
        end
    end
end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))



%% Performance Profile (Dolan & More, Math. Program. 91, 2002)

% sequential runs stored per Hessian estimation method (cf. hess_est_method)
seqMethodStructs = {seqSolverStructs_regularization};
seqMethodNames   = {hess_est_method_1};
if enable_second_sequential
    seqMethodStructs{end+1} = seqSolverStructs_hessApprox2;
    seqMethodNames{end+1}   = hess_est_method_2;
end

% solver names: one sequential solver per Hessian estimation method plus coordinate descent
% (explicit loop: arrayfun/concat can produce nested cells, which fprintf rejects)
solverNames = cell(1, numel(seqMethodNames) + 1);
for s = 1:numel(seqMethodNames)
    solverNames{s} = ['Sequential (' seqMethodNames{s} ')'];
end
solverNames{end} = 'Coordinate-descent';

% evaluate the stored test data: shifted geometric means and the
% performance profile plots (legend labels passed as inputs)
plotPerformanceProfiles(seqMethodStructs, cdSolverStructs, solverNames, commonOpts.maxSolveTime, probNames)

% remove the added folders from the path again now that the benchmarks
% are done (exactly the folders added at the top of this script)
for ii = 1:numel(addPaths)
    rmpath(addPaths{ii})
end
