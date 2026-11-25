% %% ***********************************************************************
% %
% %  Region-of-attraction Estimation Benchmark Test Runner
% %
% %  Note 1: Each runs uses various combinations of polynomial degrees. For SOS
% %  problem to be feasible or to improve solution quality higher degree are
% %  needed. Hence, with a higher degree more decision variables occur and
% %  the Gram matrices become larger. Therefore, we can check the behavior
% %  for over several benchmark problems, i.e., we have sub-benchmarks so to
% %  say that increase our benchmark set.
% %
% %  Note 2: The user might select a second hessian approximation method for
% %  sequential. This is optional, but interesting in terms of convergence
% %  speed.
% %
% %
% %
% % ***********************************************************************
% set(groot, ...
%     'defaultTextInterpreter', 'latex', ...
%     'defaultAxesTickLabelInterpreter', 'latex', ...
%     'defaultLegendInterpreter', 'latex');
% 
close all
clear
clc

% set MATLAB threads to 1
maxNumCompThreads(1);


date_of_benchmark = string(datetime("today","Format","dd-MMM-yyyy"));
name_of_benchmark = 'sensitivity_initial_guess';

% get versions
matlabVersion = version;
casadiVersion = casadi.CasadiMeta.version;


% add the problem specific folders and the shared helper functions
runnerDir = pwd;
addPaths  = { ...
    fullfile(runnerDir,'Aircraft'), ...
    fullfile(runnerDir,'../../','helperFunctions')};

for ii = 1:numel(addPaths)
    addpath(addPaths{ii})
end

% enable plotting of convergence for each run
enable_plot_convergence  = false;
enable_second_sequential = true;

% Hessian estimation methods of the sequential solver to compare
hess_est_method_1 = 'regularization_frobenius'; % our default (see also sequentialCommon.m)
hess_est_method_2 = 'damped_bfgs';      % 'damped_bfgs', 'regularization_levenberg', 'regularization_mirror', 'regularization_frobenius'

% common options shared by the methods
commonOpts.maxiter       = 500;
commonOpts.tolerance_opt = 1e-4;
commonOpts.eps_cost      = 1e-6;
commonOpts.small_prog_count = 10;
commonOpts.maxSolveTime     = 3600;
commonOpts.almost_opt_count = 10;
commonOpts.solver           =   'mosek';
commonOpts.miniter          = 5; % minimum iterations before an aborted CD run still counts as successful

% MOSEK settings

if strcmpi(commonOpts.solver,'mosek')
    commonOpts.sdp_solver_opts.mosek_param = struct('MSK_IPAR_NUM_THREADS',1);
    commonOpts.sdp_solver_opts.mosek_echo = 0;

elseif strcmpi(commonOpts.solver,'clarabel')

    commonOpts.sdp_solver_opts.direct_solve_method =  'FAER';
    commonOpts.sdp_solver_opts.verbose             = 0;

elseif strcmpi(commonOpts.solver,'scs')

    commonOpts.sdp_solver_opts.dense        = false;
    commonOpts.sdp_solver_opts.use_indirect = false;    % conjugate gradient (iterative)
    commonOpts.sdp_solver_opts.verbose      = false;

end

% variable names of the results table (one row per run)
tableVariableNames = {'Prob. No.','Problem Name','no. con. (SDP)','no. dec. (SDP)','Size Gram matrices', ...
    'Build Time (CD)','Build Time (SQ)','Solve Time (CD)','Solve Time (SQ)', ...
    'Iterations (CD)','Iterations (SQ)','time/iter. (CD)','time/iter. (SQ)','Final Cost (CD)','Final Cost (SQ)', ...
    'Solved (CD)','Solved (SQ)'};

% extend the results table by the sequential run with the second Hessian
% estimation method (cf. enable_second_sequential)
if enable_second_sequential
    tableVariableNames = [tableVariableNames, ...
        ['Build Time (SQ ' hess_est_method_2 ')'], ...
        ['Solve Time (SQ ' hess_est_method_2 ')'], ...
        ['Iterations (SQ ' hess_est_method_2 ')'], ...
        ['time/iter. (SQ ' hess_est_method_2 ')'], ...
        ['Final Cost (SQ ' hess_est_method_2 ')'], ...
        ['Solved (SQ ' hess_est_method_2 ')']];
end

% initialize the cell array that accumulates one row per run of the
% results table (one column per table variable)
tableRows = cell(0, numel(tableVariableNames));


% initialize the cell arrays of stored runs (one entry per run)
seqSolverStructs_regularization = {};
seqSolverStructs_hessApprox2    = {};
cdSolverStructs                 = {};
probNames                       = {}; % problem name per run (for the evaluation output)


total_benchmark_time_start = tic;

% 
%% Aircraft 4D

% define the problem instances based on different polynomial degrees
degV = {2,   2, 2:4, 1:4};
degS = {2, 2:4, 2:4, 2:4};


degV = {2};
degS = {2:4};

for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Aircraft = generateAircraft_probStruct(degV{kk},degS{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Initial Guess Sensitivity for %s Instance %d of %d. \n', probStruct_Aircraft.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_Aircraft(probStruct_Aircraft);
    
    for kkk = 1:numel(x0)
    % solve Aircraft example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Aircraft_regularization = solve_ROA_seq_bench_fun(probStruct_Aircraft, ...
                                                                    commonOpts, ...
                                                                    x0{kkk}, ...
                                                                    hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_Aircraft_regularization;

    % extract size and timing statistic from sequential
    evalSequentialStruct_Aircraft_regularization                 = evalSequential(solStruct_seq_Aircraft_regularization);
    evalSequentialStruct_Aircraft_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Aircraft_regularization.probNumber      = kk*kkk;

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_Aircraft_hessApprox2 = solve_ROA_seq_bench_fun(probStruct_Aircraft, ...
                                                                     commonOpts, ...
                                                                     x0{kkk}, ...
                                                                     hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_Aircraft_hessApprox2;

        % just some gap in display output
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        % extract size and timing statistic from sequential
        evalSequentialStruct_Aircraft_hessApprox2                 = evalSequential(solStruct_seq_Aircraft_hessApprox2);
        evalSequentialStruct_Aircraft_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_Aircraft_hessApprox2.probNumber      = kk*kkk;
    end

    % solve Aircraft example with two states using custom coordinate-descent
    solStruct_CD_Aircraft = solve_ROA_CD_bench_fun(probStruct_Aircraft, ...
                                                   commonOpts, ...
                                                   x0{kkk});

    cdSolverStructs{end+1} = solStruct_CD_Aircraft;
    probNames{end+1}       = sprintf('%s (instance %d of %d)', probStruct_Aircraft.Name, kk, numel(degV));

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % Generate table

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows, ...
                                    probStruct_Aircraft, ...
                                    solStruct_CD_Aircraft, ...
                                    solStruct_seq_Aircraft_regularization, ...
                                    solStruct_seq_Aircraft_hessApprox2, ...
                                    true);
    else
        tableRows = appendResultRow(tableRows, ...
                                    probStruct_Aircraft, ...
                                    solStruct_CD_Aircraft, ...
                                    solStruct_seq_Aircraft_regularization, ...
                                    [], ...
                                    true);
    end

    if enable_plot_convergence
        if enable_second_sequential
            % Evaluate
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_Aircraft_regularization,...
                            evalSequentialStruct_Aircraft_hessApprox2,...
                            solStruct_CD_Aircraft)
        else
            % Evaluate
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_Aircraft_regularization,...
                            [],...
                            solStruct_CD_Aircraft)
        end
    end
    end
end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))




total_benchmark_time = toc(total_benchmark_time_start);
total_benchmark_time_h = total_benchmark_time/3600;

fprintf('Total benchmark time: %f second\n', total_benchmark_time);
fprintf('Total benchmark time: %f minutes\n', total_benchmark_time/60);
fprintf('Total benchmark time: %f hours\n', total_benchmark_time_h);


% %% Store benchmark data
% 
% % fprintf('\n')
% % fprintf('Store benchmark data in .mat file\n')
% % fprintf('\n')
% % 
% % save([name_of_benchmark '_' commonOpts.solver '.mat'],...         % name of .mat files
% %      "name_of_benchmark",... 
% %          "total_benchmark_time_h",...
% %      'casadiVersion',...                    % store the used casadi version
% %      'date_of_benchmark',...                % date when the benchmark was generated
% %      'matlabVersion',...                    % store the matlab version
% %      'probNames',...                        % names for each run
% %      'cdSolverStructs',...                  % struct containing stats and solution for CD
% %      'seqSolverStructs_hessApprox2',...     % struct containing stats and solution for the sequential with alternative hessian estimation
% %      'seqSolverStructs_regularization',...  % struct containing stats and solution for the sequential with default hessian estimation
% %      'commonOpts',...                       % shared options
% %      'enable_second_sequential', ...        % flag to enable second sequential run with alternative hessian approx
% %      'enable_plot_convergence', ...         % flag to enable plot of convergence
% %      'hess_est_method_1', ...               % default hessian estimation method (needed as string for plotting)
% %      'hess_est_method_2', ...               % second hessian estimation method
% %      'tableRows', ...                       % table with result
% %      'tableVariableNames')                  % column names for table
% % 
% % fprintf('\n')
% 
 %% Performance Profile (Dolan & More, Math. Program. 91, 2002)


fprintf('Evaluate benchmark ... \n')
fprintf('\n')
% sequential runs stored per Hessian estimation method (cf. hess_est_method)
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

% evaluate the stored test data: shifted geometric means and the
% performance profile plots (legend labels passed as inputs)
plotPerformanceProfiles(seqMethodStructs, ...
                        cdSolverStructs, ...
                        solverNames, ...
                        commonOpts.maxSolveTime, ...
                        probNames)


% final results table accumulated over all runs and all problems
resultsTable = cell2table(tableRows,'VariableNames',tableVariableNames);
disp(resultsTable)


% remove the added folders from the path again now that the benchmarks
% % are done (exactly the folders added at the top of this script)
% for ii = 1:numel(addPaths)
%     rmpath(addPaths{ii})
% end


%% Eval cost
fprintf('Evaluate benchmark ... \n')
fprintf('\n')
% sequential runs stored per Hessian estimation method (cf. hess_est_method)
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

cost_spread_calculation(seqMethodStructs, ...
    cdSolverStructs, ...
    solverNames, ...
    commonOpts.maxSolveTime, ...
    probNames)


% get minimum cost overall feasible runs for each method

% compute spread
