%% ***********************************************************************
%
%  Reachability Estimation Benchmark Runner
%
%
%   SOS Problem for Reachability Estimation
%
%   min      ||g - (V(0)-β)||
%   V,s,κ,β 
%
%   s.j.    s_i                                 is SOS
%           β                                    is SOS
%           s1(V-β) - s2 h_T - ∂V/∂t - ∇Vf(x,κ)  is SOS
%           s3(V-β) - s4 h_T + κ - u_min         is SOS
%           s5(V-β) - s6 h_T + u_max - κ         is SOS
%           s7(V-β) - s8 h_T - g                 is SOS
%           s9(V(T)-β) - g_T                     is SOS
%
%   V    - Reachability function
%   s_i  - SOS multipliers
%   κ    - Control law
%   β    - Level set
%   h_T  - Time polynomial t(T-t)
%   g    - Safe set
%   g_T  - Terminal set
%   T    - Time horizon
%   u_min, u_max - Control bounds
%  
%
%  Overview
%  
%  Benchmark 1: VDP
%  Benchmark 2: MIMO
%  Benchmark 3: Aircraft
%  Benchmark 4: Satellite
%
%
%  License: GNU GENERAL PUBLIC LICENSE Version 3
%
%% ***********************************************************************

set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');

close all
clear
clc


% set MATLAB threads to 1
maxNumCompThreads(1);

date_of_benchmark = string(datetime("today","Format","dd-MMM-yyyy"));
name_of_benchmark = 'reachability_benchmark';

% get versions
matlabVersion = version;
casadiVersion = casadi.CasadiMeta.version;

% add the problem specific folders and the shared helper functions
runnerDir = pwd;

addPaths  = { ...
    fullfile(runnerDir,'vdp_reachability'), ...
    fullfile(runnerDir,'mimo_reachability'), ...
    fullfile(runnerDir,'aircraft_reachability'), ...
    fullfile(runnerDir,'satellite_reachability'),...
    fullfile(runnerDir,'..','helperFunctions')};

for ii = 1:numel(addPaths)
    addpath(addPaths{ii})
end

% enable plotting of convergence for each run
enable_plot_convergence  = true;
enable_second_sequential = true;

% Hessian estimation methods of the sequential solver to compare
hess_est_method_1 = 'regularization_frobenius'; % our default (see also SequentialCommon.m)
hess_est_method_2 = 'damped_bfgs';              % 'damped_bfgs', 'regularization_levenberg', 'regularization_mirror', 'regularization_frobenius'

% common options shared by the methods
commonOpts.maxiter          = 500;
commonOpts.tolerance_opt    = 1e-4;
commonOpts.eps_cost         = 1e-6;
commonOpts.small_prog_count = 10;
commonOpts.maxSolveTime     = 3600;
commonOpts.almost_opt_count = 10;
commonOpts.solver           = 'mosek';
commonOpts.miniter          = 5; % minimum iterations before an aborted CD run still counts as successful

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
degV = {0:2 , 2:4, 0:4};
degK = {0:1 , 2:4, 0:4};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_VDP = generate_prob_struct_vdp_reachability(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Reachbility Analysis Benchmark for %s Instance %d of %d. \n', probStruct_VDP.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_vdp_reachability(probStruct_VDP);

    % solve VDP example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_VDP_regularization = solve_reachability_seq_bench_fun(probStruct_VDP, ...
                                                                        commonOpts, ...
                                                                        x0, ...
                                                                        hess_est_method);

    % store in struct
    seqSolverStructs_regularization{end+1} = solStruct_seq_VDP_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_VDP_regularization                 = evalSequential(solStruct_seq_VDP_regularization);
    evalSequentialStruct_VDP_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_VDP_regularization.probNumber      = kk;

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % if enabled compute with second hessian estimation method
    if enable_second_sequential

        solStruct_seq_VDP_hessApprox2 = solve_reachability_seq_bench_fun(probStruct_VDP, ...
                                                                         commonOpts, ...
                                                                         x0, ...
                                                                         hess_est_method_2);

        % store in struct
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
    else
        evalSequentialStruct_VDP_hessApprox2 = [];
        solStruct_seq_VDP_hessApprox2        = [];
    end

    % solve VDP example with two states using custom coordinate-descent
    solStruct_CD_VDP = solve_reachability_CD_bench_fun(probStruct_VDP, ...
                                                       commonOpts, ...
                                                       x0);

    cdSolverStructs{end+1} = solStruct_CD_VDP;
    probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_VDP.Name, kk, numel(degV), commonOpts.solver);

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows, ...
                                probStruct_VDP, ...
                                solStruct_CD_VDP, ...
                                solStruct_seq_VDP_regularization, ...
                                solStruct_seq_VDP_hessApprox2, ...
                                true);

    % plot the scaled first order KKT condition over iterations
    if enable_plot_convergence
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_VDP_regularization,...
                            evalSequentialStruct_VDP_hessApprox2,...
                            solStruct_CD_VDP)
    end
end % end of for-loop for VDP benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))

%% MIMO

% define the problem instances based on different polynomial degrees
degV = {0:2 , 2:4, 0:4};
degK = {0:1 , 2:4, 0:4};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_MIMO = generate_prob_struct_mimo_reachability(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Reachbility Analysis Benchmark for %s Instance %d of %d. \n', probStruct_MIMO.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_mimo_reachability(probStruct_MIMO);

    % solve MIMO example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_MIMO_regularization = solve_reachability_seq_bench_fun(probStruct_MIMO, ...
                                                                        commonOpts, ...
                                                                        x0, ...
                                                                        hess_est_method);

    % store in struct
    seqSolverStructs_regularization{end+1} = solStruct_seq_MIMO_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_MIMO_regularization                 = evalSequential(solStruct_seq_MIMO_regularization);
    evalSequentialStruct_MIMO_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_MIMO_regularization.probNumber      = kk;

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % if enabled compute with second hessian estimation method
    if enable_second_sequential

        solStruct_seq_MIMO_hessApprox2 = solve_reachability_seq_bench_fun(probStruct_MIMO, ...
                                                                         commonOpts, ...
                                                                         x0, ...
                                                                         hess_est_method_2);

        % store in struct
        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_MIMO_hessApprox2;

        % just some gap in display output
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        % extract size and timing statistic from sequential
        evalSequentialStruct_MIMO_hessApprox2                 = evalSequential(solStruct_seq_MIMO_hessApprox2);
        evalSequentialStruct_MIMO_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_MIMO_hessApprox2.probNumber      = kk;
    else
        evalSequentialStruct_MIMO_hessApprox2 = [];
        solStruct_seq_MIMO_hessApprox2        = [];
    end

    % solve MIMO example with two states using custom coordinate-descent
    solStruct_CD_MIMO = solve_reachability_CD_bench_fun(probStruct_MIMO, ...
                                                        commonOpts, ...
                                                        x0);

    cdSolverStructs{end+1} = solStruct_CD_MIMO;
    probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_MIMO.Name, kk, numel(degV), commonOpts.solver);

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows, ...
                                probStruct_MIMO, ...
                                solStruct_CD_MIMO, ...
                                solStruct_seq_MIMO_regularization, ...
                                solStruct_seq_MIMO_hessApprox2, ...
                                true);

    % plot the scaled first order KKT condition over iterations
    if enable_plot_convergence
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_MIMO_regularization,...
                            evalSequentialStruct_MIMO_hessApprox2,...
                            solStruct_CD_MIMO)
    end
end % end of for-loop for MIMO benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))

%% Aircraft

degV = {1:2 , 2:4};
degK = {0:1 , 0:1};


for kk = 1:numel(degV)

    probStruct_Aircraft = generate_prob_struct_aircraft_reachability(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Reachbility Analysis Benchmark for %s Instance %d of %d. \n', probStruct_Aircraft.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    x0 = generate_init_guess_aircraft_reachability(probStruct_Aircraft);

    hess_est_method = hess_est_method_1;
    solStruct_seq_Aircraft_regularization = solve_reachability_seq_bench_fun(probStruct_Aircraft, ...
                                                                        commonOpts, ...
                                                                        x0, ...
                                                                        hess_est_method);

    % store in struct
    seqSolverStructs_regularization{end+1} = solStruct_seq_Aircraft_regularization;

    evalSequentialStruct_Aircraft_regularization                 = evalSequential(solStruct_seq_Aircraft_regularization);
    evalSequentialStruct_Aircraft_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Aircraft_regularization.probNumber      = kk;

    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % if enabled compute with second hessian estimation method
    if enable_second_sequential

        solStruct_seq_Aircraft_hessApprox2 = solve_reachability_seq_bench_fun(probStruct_Aircraft, ...
                                                                         commonOpts, ...
                                                                         x0, ...
                                                                         hess_est_method_2);

        % store in struct
        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_Aircraft_hessApprox2;

        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        evalSequentialStruct_Aircraft_hessApprox2                 = evalSequential(solStruct_seq_Aircraft_hessApprox2);
        evalSequentialStruct_Aircraft_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_Aircraft_hessApprox2.probNumber      = kk;
    else
        evalSequentialStruct_Aircraft_hessApprox2 = [];
        solStruct_seq_Aircraft_hessApprox2        = [];
    end

    solStruct_CD_Aircraft = solve_reachability_CD_bench_fun(probStruct_Aircraft, ...
                                                            commonOpts, ...
                                                            x0);

    cdSolverStructs{end+1} = solStruct_CD_Aircraft;
    probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_Aircraft.Name, kk, numel(degV), commonOpts.solver);

    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    tableRows = appendResultRow(tableRows, ...
                                probStruct_Aircraft, ...
                                solStruct_CD_Aircraft, ...
                                solStruct_seq_Aircraft_regularization, ...
                                solStruct_seq_Aircraft_hessApprox2, ...
                                true);

    % plot the scaled first order KKT condition over iterations
    if enable_plot_convergence
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_Aircraft_regularization,...
                            evalSequentialStruct_Aircraft_hessApprox2,...
                            solStruct_CD_Aircraft)
    end
end % end of for-loop for Aircraft benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Satellite

degV = {1:2 , 2:4};
degK = {0:1 , 0:1};


for kk = 1:numel(degV)

    probStruct_Satellite = generate_prob_struct_satellite_reachability(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Reachbility Analysis Benchmark for %s Instance %d of %d. \n', probStruct_Satellite.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    x0 = generate_init_guess_satellite_reachability(probStruct_Satellite);

    hess_est_method = hess_est_method_1;
    solStruct_seq_Satellite_regularization = solve_reachability_seq_bench_fun(probStruct_Satellite, ...
        commonOpts, ...
        x0, ...
        hess_est_method);

    % store in struct
    seqSolverStructs_regularization{end+1} = solStruct_seq_Satellite_regularization;

    evalSequentialStruct_Satellite_regularization                 = evalSequential(solStruct_seq_Satellite_regularization);
    evalSequentialStruct_Satellite_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Satellite_regularization.probNumber      = kk;

    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % if enabled compute with second hessian estimation method
    if enable_second_sequential

        solStruct_seq_Satellite_hessApprox2 = solve_reachability_seq_bench_fun(probStruct_Satellite, ...
            commonOpts, ...
            x0, ...
            hess_est_method_2);

        % store in struct
        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_Satellite_hessApprox2;

        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        evalSequentialStruct_Satellite_hessApprox2                 = evalSequential(solStruct_seq_Satellite_hessApprox2);
        evalSequentialStruct_Satellite_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_Satellite_hessApprox2.probNumber      = kk;
    else
        evalSequentialStruct_Satellite_hessApprox2 = [];
        solStruct_seq_Satellite_hessApprox2        = [];
    end

    solStruct_CD_Satellite = solve_reachability_CD_bench_fun(probStruct_Satellite, ...
                                                             commonOpts, ...
                                                             x0);

    cdSolverStructs{end+1} = solStruct_CD_Satellite;
    probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_Satellite.Name, kk, numel(degV), commonOpts.solver);

    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    tableRows = appendResultRow(tableRows, ...
        probStruct_Satellite, ...
        solStruct_CD_Satellite, ...
        solStruct_seq_Satellite_regularization, ...
        solStruct_seq_Satellite_hessApprox2, ...
        true);

    % plot the scaled first order KKT condition over iterations
    if enable_plot_convergence
            plotConvergence(probNames{end}, ...
                enable_second_sequential,...
                evalSequentialStruct_Satellite_regularization,...
                evalSequentialStruct_Satellite_hessApprox2,...
                solStruct_CD_Satellite)
    end
end % end of for-loop for Satellite benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))



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

% solver names: one sequential solver per Hessian estimation method plus coordinate descent
% (explicit loop: arrayfun/concat can produce nested cells, which fprintf rejects)
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


% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))



% remove the added folders from the path again now that the benchmarks
% are done (exactly the folders added at the top of this script)
for ii = 1:numel(addPaths)
    rmpath(addPaths{ii})
end
