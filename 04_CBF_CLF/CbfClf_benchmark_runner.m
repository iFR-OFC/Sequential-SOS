%% ***********************************************************************
%
%  CBF-CLF Benchmark Test Runner
% ***********************************************************************
set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');

close all
clear
clc


maxNumCompThreads(1);

date_of_benchmark = string(datetime("today","Format","dd-MMM-yyyy"));
name_of_benchmark = 'CBF_CLF_benchmark';

% get versions
matlabVersion = version;
casadiVersion = casadi.CasadiMeta.version;

% add the problem specific folders and the shared helper functions
runnerDir = pwd;

addPaths  = { ...
    fullfile(runnerDir,'Satellite'), ...
    fullfile(runnerDir,'VDP/'), ...
    fullfile(runnerDir,'VDP_coupled/'), ...
   fullfile(runnerDir,'Aircraft_affine/'), ...
    fullfile(runnerDir,'..','helperFunctions')};

for ii = 1:numel(addPaths)
    addpath(addPaths{ii})
end

% enable plotting of convergence for each run
enable_plot_convergence  = true;
enable_second_sequential = true;

% Hessian estimation methods of the sequential solver to compare
hess_est_method_1 = 'regularization_frobenius'; % our default (see also sequentialCommon.m)
hess_est_method_2 = 'regularization_levenberg';              % 'damped_bfgs', 'regularization_levenberg', 'regularization_mirror', 'regularization_frobenius'

% common options shared by the methods
commonOpts.maxiter       = 500;
commonOpts.tolerance_opt = 1e-4;
commonOpts.eps_cost      = 1e-6;
commonOpts.small_prog_count = 10;
commonOpts.maxSolveTime  = 3600;
commonOpts.almost_opt_count = 10;
commonOpts.solver        = 'mosek';
commonOpts.miniter       = 5; % minimum iterations before an aborted CD run still counts as successful

% MOSEK settings

if strcmpi(commonOpts.solver,'mosek')
    commonOpts.sdp_solver_opts.mosek_param = struct( 'MSK_IPAR_NUM_THREADS', 1);
    commonOpts.sdp_solver_opts.mosek_echo  = 0;

elseif strcmpi(commonOpts.solver,'clarabel')
    commonOpts.sdp_solver_opts.direct_solve_method =  'FAER';
    commonOpts.sdp_solver_opts.verbose             = 0;
elseif strcmpi(commonOpts.solver,'scs')
    commonOpts.sdp_solver_opts.dense               = 1;
    commonOpts.sdp_solver_opts.verbose             = 0;
end


% variable names of the results table (one row per run)
tableVariableNames = {'Prob. No.','Problem Name','no. con. (SDP)','no. dec. (SDP)','Size Gram matrices', ...
                        'Build Time (CD)','Build Time (SQ)','Solve Time (CD)','Solve Time (SQ)', ...
                        'Iterations (CD)','Iterations (SQ)','time/iter. (CD)','time/iter. (SQ)','Final Cost (CD)','Final Cost (SQ)', ...
                        'Solved (CD)','Solved (SQ)'};

% extend the results table by the sequential run with the second Hessian
% estimation method (if enable_second_sequential == 1)
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
% results table
tableRows = cell(0, numel(tableVariableNames));

% initialize the cell arrays of stored runs (one entry per run)
seqSolverStructs_regularization = {};
seqSolverStructs_hessApprox2    = {};
cdSolverStructs                 = {};
probNames                       = {}; % problem name per run (for the evaluation output)

total_benchmark_time_start = tic;


%% VDP

% define the problem instances based on different polynomial degrees
degV = {2,   2, 2, 2:4};
degH = {2, 2:4, 2, 2:4};
degK = {1, 1, 1:2,  1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_VDP = generateVDP_probStruct(degV{kk},degH{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('CBF-CLF Benchmark for %s Instance %d of %d. \n', probStruct_VDP.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_VDP(probStruct_VDP);

    % solve VDP example using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_VDP_regularization = solve_CbfClf_seq_bench_fun(probStruct_VDP, ...
                                                                        commonOpts, ...
                                                                        x0, ...
                                                                        hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_VDP_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_VDP_regularization                 = evalSequential(solStruct_seq_VDP_regularization);
    evalSequentialStruct_VDP_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_VDP_regularization.probNumber      = kk;
    evalSequentialStruct_VDP_hessApprox2 = [];

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_VDP_hessApprox2 = solve_CbfClf_seq_bench_fun(probStruct_VDP, ...
                                                                         commonOpts, ...
                                                                         x0, ...
                                                                         hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_VDP_hessApprox2;

        % just some gap in display output
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        % extract size and timing statistic from sequential
        evalSequentialStruct_VDP_hessApprox2                 = evalSequential(solStruct_seq_VDP_hessApprox2);
        evalSequentialStruct_VDP_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_VDP_hessApprox2.probNumber      = kk;
    end

    % solve VDP example using custom coordinate-descent
    solStruct_CD_VDP = solve_CbfClf_CD_bench_fun(probStruct_VDP,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_VDP;
    probNames{end+1}       = sprintf('%s (instance %d of %d)', probStruct_VDP.Name, kk, numel(degV));

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows, ...
                                    probStruct_VDP, ...
                                    solStruct_CD_VDP, ...
                                    solStruct_seq_VDP_regularization, ...
                                    solStruct_seq_VDP_hessApprox2, ...
                                    true);
    else
        tableRows = appendResultRow(tableRows, ...
                                    probStruct_VDP, ...
                                    solStruct_CD_VDP, ...
                                    solStruct_seq_VDP_regularization, ...
                                    [], ...
                                    true);
    end

    if enable_plot_convergence
        if enable_second_sequential
            % Evaluate
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_VDP_regularization,...
                            evalSequentialStruct_VDP_hessApprox2,...
                            solStruct_CD_VDP)
        else
            % Evaluate
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_VDP_regularization,...
                            [],...
                            solStruct_CD_VDP)
        end
    end
end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))

%% VDP_coupled

% define the problem instances based on different polynomial degrees
degV = {2,   2, 2, 2:4};
degH = {2, 2:4, 2, 2:4};
degK = {1, 1, 1:2,  1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_VDP_coupled = generateVDP_coupled_probStruct(degV{kk},degH{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('CBF-CLF Benchmark for %s Instance %d of %d. \n', probStruct_VDP_coupled.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_VDP_coupled(probStruct_VDP_coupled);

    % solve VDP_coupled example using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_VDP_coupled_regularization = solve_CbfClf_seq_bench_fun(probStruct_VDP_coupled, ...
                                                                        commonOpts, ...
                                                                        x0, ...
                                                                        hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_VDP_coupled_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_VDP_coupled_regularization                 = evalSequential(solStruct_seq_VDP_coupled_regularization);
    evalSequentialStruct_VDP_coupled_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_VDP_coupled_regularization.probNumber      = kk;
    evalSequentialStruct_VDP_coupled_hessApprox2 = [];

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_VDP_coupled_hessApprox2 = solve_CbfClf_seq_bench_fun(probStruct_VDP_coupled, ...
                                                                         commonOpts, ...
                                                                         x0, ...
                                                                         hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_VDP_coupled_hessApprox2;

        % just some gap in display output
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        % extract size and timing statistic from sequential
        evalSequentialStruct_VDP_coupled_hessApprox2                 = evalSequential(solStruct_seq_VDP_coupled_hessApprox2);
        evalSequentialStruct_VDP_coupled_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_VDP_coupled_hessApprox2.probNumber      = kk;
    end

    % solve VDP_coupled example using custom coordinate-descent
    solStruct_CD_VDP_coupled = solve_CbfClf_CD_bench_fun(probStruct_VDP_coupled,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_VDP_coupled;
    probNames{end+1}       = sprintf('%s (instance %d of %d)', probStruct_VDP_coupled.Name, kk, numel(degV));

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows, ...
                                    probStruct_VDP_coupled, ...
                                    solStruct_CD_VDP_coupled, ...
                                    solStruct_seq_VDP_coupled_regularization, ...
                                    solStruct_seq_VDP_coupled_hessApprox2, ...
                                    true);
    else
        tableRows = appendResultRow(tableRows, ...
                                    probStruct_VDP_coupled, ...
                                    solStruct_CD_VDP_coupled, ...
                                    solStruct_seq_VDP_coupled_regularization, ...
                                    [], ...
                                    true);
    end

    if enable_plot_convergence
        if enable_second_sequential
            % Evaluate
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_VDP_coupled_regularization,...
                            evalSequentialStruct_VDP_coupled_hessApprox2,...
                            solStruct_CD_VDP_coupled)
        else
            % Evaluate
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_VDP_coupled_regularization,...
                            [],...
                            solStruct_CD_VDP_coupled)
        end
    end
end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Aircraft_affine

% define the problem instances based on different polynomial degrees
degV = {2,   2, 2, 2:4};
degH = {2, 2:4, 2, 2:4};
degK = {1, 1, 1:2,  1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Aircraft_affine = generateAircraft_affine_probStruct(degV{kk},degH{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('CBF-CLF Benchmark for %s Instance %d of %d. \n', probStruct_Aircraft_affine.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_Aircraft_affine(probStruct_Aircraft_affine);

    % solve Aircraft_affine example using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Aircraft_affine_regularization = solve_CbfClf_seq_bench_fun(probStruct_Aircraft_affine, ...
                                                                        commonOpts, ...
                                                                        x0, ...
                                                                        hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_Aircraft_affine_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_Aircraft_affine_regularization                 = evalSequential(solStruct_seq_Aircraft_affine_regularization);
    evalSequentialStruct_Aircraft_affine_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Aircraft_affine_regularization.probNumber      = kk;
    evalSequentialStruct_Aircraft_affine_hessApprox2 = [];

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_Aircraft_affine_hessApprox2 = solve_CbfClf_seq_bench_fun(probStruct_Aircraft_affine, ...
                                                                         commonOpts, ...
                                                                         x0, ...
                                                                         hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_Aircraft_affine_hessApprox2;

        % just some gap in display output
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        % extract size and timing statistic from sequential
        evalSequentialStruct_Aircraft_affine_hessApprox2                 = evalSequential(solStruct_seq_Aircraft_affine_hessApprox2);
        evalSequentialStruct_Aircraft_affine_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_Aircraft_affine_hessApprox2.probNumber      = kk;
    end

    % solve Aircraft_affine example using custom coordinate-descent
    solStruct_CD_Aircraft_affine = solve_CbfClf_CD_bench_fun(probStruct_Aircraft_affine,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_Aircraft_affine;
    probNames{end+1}       = sprintf('%s (instance %d of %d)', probStruct_Aircraft_affine.Name, kk, numel(degV));

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows, ...
                                    probStruct_Aircraft_affine, ...
                                    solStruct_CD_Aircraft_affine, ...
                                    solStruct_seq_Aircraft_affine_regularization, ...
                                    solStruct_seq_Aircraft_affine_hessApprox2, ...
                                    true);
    else
        tableRows = appendResultRow(tableRows, ...
                                    probStruct_Aircraft_affine, ...
                                    solStruct_CD_Aircraft_affine, ...
                                    solStruct_seq_Aircraft_affine_regularization, ...
                                    [], ...
                                    true);
    end

    if enable_plot_convergence
        if enable_second_sequential
            % Evaluate
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_Aircraft_affine_regularization,...
                            evalSequentialStruct_Aircraft_affine_hessApprox2,...
                            solStruct_CD_Aircraft_affine)
        else
            % Evaluate
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_Aircraft_affine_regularization,...
                            [],...
                            solStruct_CD_Aircraft_affine)
        end
    end
end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Satellite

% define the problem instances based on different polynomial degrees
degV = {2,2,2:4};
degH = {2,2:4,2};
degK = {1,1:2,1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Satellite = generateSatellite_probStruct(degV{kk},degH{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('CBF-CLF Benchmark for %s Instance %d of %d. \n', probStruct_Satellite.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_Satellite(probStruct_Satellite);

    % solve Satellite example using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Satellite_regularization = solve_CbfClf_seq_bench_fun(probStruct_Satellite, ...
                                                                        commonOpts, ...
                                                                        x0, ...
                                                                        hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_Satellite_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_Satellite_regularization                 = evalSequential(solStruct_seq_Satellite_regularization);
    evalSequentialStruct_Satellite_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Satellite_regularization.probNumber      = kk;
    evalSequentialStruct_Satellite_hessApprox2 = [];

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_Satellite_hessApprox2 = solve_CbfClf_seq_bench_fun(probStruct_Satellite, ...
                                                                         commonOpts, ...
                                                                         x0, ...
                                                                         hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_Satellite_hessApprox2;

        % just some gap in display output
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        % extract size and timing statistic from sequential
        evalSequentialStruct_Satellite_hessApprox2                 = evalSequential(solStruct_seq_Satellite_hessApprox2);
        evalSequentialStruct_Satellite_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_Satellite_hessApprox2.probNumber      = kk;
    end

    % solve Satellite example using custom coordinate-descent
    solStruct_CD_Satellite = solve_CbfClf_CD_bench_fun(probStruct_Satellite,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_Satellite;
    probNames{end+1}       = sprintf('%s (instance %d of %d)', probStruct_Satellite.Name, kk, numel(degV));

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows, ...
                                    probStruct_Satellite, ...
                                    solStruct_CD_Satellite, ...
                                    solStruct_seq_Satellite_regularization, ...
                                    solStruct_seq_Satellite_hessApprox2, ...
                                    true);
    else
        tableRows = appendResultRow(tableRows, ...
                                    probStruct_Satellite, ...
                                    solStruct_CD_Satellite, ...
                                    solStruct_seq_Satellite_regularization, ...
                                    [], ...
                                    true);
    end

    if enable_plot_convergence
        if enable_second_sequential
            % Evaluate
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_Satellite_regularization,...
                            evalSequentialStruct_Satellite_hessApprox2,...
                            solStruct_CD_Satellite)
        else
            % Evaluate
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_Satellite_regularization,...
                            [],...
                            solStruct_CD_Satellite)
        end
    end
end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Results table (all benchmarks, all runs)

% final results table accumulated over all runs and all problems
resultsTable = cell2table(tableRows,'VariableNames',tableVariableNames);
disp(resultsTable)


total_benchmark_time = toc(total_benchmark_time_start);
total_benchmark_time_h = total_benchmark_time/3600;

fprintf('Total benchmark time: %f second\n', total_benchmark_time);
fprintf('Total benchmark time: %f minutes\n', total_benchmark_time/60);
fprintf('Total benchmark time: %f hours\n', total_benchmark_time_h);

%% Store benchmark data

fprintf('\n')
fprintf('Store benchmark data in .mat file\n')
fprintf('\n')

save([name_of_benchmark '_' commonOpts.solver '.mat'],...         % name of .mat files
         "name_of_benchmark",... 
         "total_benchmark_time_h",...
     'casadiVersion',...                    % store the used casadi version
     'date_of_benchmark',...                % date when the benchmark was generated
     'matlabVersion',...                    % store the matlab version
     'probNames',...                        % names for each run
     'cdSolverStructs',...                  % struct containing stats and solution for CD
     'seqSolverStructs_hessApprox2',...     % struct containing stats and solution for the sequential with alternative hessian estimation
     'seqSolverStructs_regularization',...  % struct containing stats and solution for the sequential with default hessian estimation
     'commonOpts',...                       % shared options
     'enable_second_sequential', ...        % flag to enable second sequential run with alternative hessian approx
     'enable_plot_convergence', ...         % flag to enable plot of convergence
     'hess_est_method_1', ...               % default hessian estimation method (needed as string for plotting)
     'hess_est_method_2', ...               % second hessian estimation method
     'tableRows', ...                       % table with result
     'tableVariableNames')                  % column names for table

fprintf('\n')

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

% solver names: one sequential solver per Hessian estimation method plus
% coordinate descent (explicit loop: arrayfun/concat can produce nested
% cells, which fprintf rejects)
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

% remove the added folders from the path again now that the benchmarks
% are done (exactly the folders added at the top of this script)
for ii = 1:numel(addPaths)
    rmpath(addPaths{ii})
end
