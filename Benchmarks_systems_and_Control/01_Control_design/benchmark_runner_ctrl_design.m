%% ***********************************************************************
%
%  Control Design Benchmark Runner
%
%
%   SOS Problem for Control Design
%
%   min      ||g - (V-β)||
%   V,s,κ,β 
%
%   s.j.    s_i                       is SOS
%           V  - ϵ(x)                 is SOS
%           β  - ϵ                    is SOS
%           s1(V-β) - ∇Vf(x,κ) - ϵ(x) is SOS
%           s2(V-β) + κ - u_min       is SOS
%           s3(V-β) + u_max - κ       is SOS
%           s4(V-β) - g               is SOS
%
%   V    - Lyapunov function
%   s_i  - SOS multipliers
%   κ    - Control law
%   β    - Level set
%   ϵ(x) - Small term
%   g    - Safe set/desired operating region
%   u_min, u_max - Control bounds
%  
%
%  Overview
%  
%  Benchmark 1: VDP
%  Benchmark 2: Lotka
%  Benchmark 3: JetEngine
%  Benchmark 4: MIMO
%  Benchmark 5: DoublePendulum
%  Benchmark 6: Aircraft (affine)
%  Benchmark 7: Pendubot
%  Benchmark 8: Spacecraft_orbit
%  Benchmark 9: Satellite
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
name_of_benchmark = 'controlDesign_benchmark';

% get versions
matlabVersion = version;
casadiVersion = casadi.CasadiMeta.version;

% add the problem specific folders and the shared helper functions
runnerDir = pwd;

addPaths  = { ...
    fullfile(runnerDir,'vdp_ctrl_design'), ...
    fullfile(runnerDir,'mimo_ctrl_design'), ...
    fullfile(runnerDir,'aircraft_affine_ctrl_design'), ...
    fullfile(runnerDir,'pendubot_ctrl_design'),...
    fullfile(runnerDir,'lotka_voltera_ctrl_design'),...
    fullfile(runnerDir,'satellite_ctrl_design'), ...
    fullfile(runnerDir,'../..','helperFunctions')};

for ii = 1:numel(addPaths)
    addpath(addPaths{ii})
end

% enable plotting of convergence for each run
enable_plot_convergence  = true; 

% enable second run of sequential SOS with different Hessian estimation
enable_second_sequential = true;

% Hessian estimation methods of the sequential solver to compare
hess_est_method_1 = 'regularization_frobenius'; % our default (see also sequentialCommon.m)
hess_est_method_2 =  'damped_bfgs';          % 'damped_bfgs', 'regularization_levenberg', 'regularization_mirror', 'regularization_frobenius'

% common options shared by the methods
commonOpts.maxiter       = 500;
commonOpts.tolerance_opt = 1e-4;
commonOpts.eps_cost      = 1e-6;
commonOpts.small_prog_count = 10;
commonOpts.maxSolveTime  =  3600;
commonOpts.almost_opt_count = 10;
commonOpts.solver        = 'mosek';
commonOpts.miniter       = 5; % minimum iterations before an aborted CD run still counts as successful

% Solver settings
if strcmpi(commonOpts.solver,'mosek')

    commonOpts.sdp_solver_opts.mosek_param = struct('MSK_IPAR_NUM_THREADS', 1); 
    commonOpts.sdp_solver_opts.mosek_echo = 0;

elseif strcmpi(commonOpts.solver,'clarabel')

    commonOpts.sdp_solver_opts.direct_solve_method =  'FAER';
    commonOpts.sdp_solver_opts.verbose             = 0;

elseif strcmpi(commonOpts.solver,'scs')

    commonOpts.sdp_solver_opts.dense               = 0;
    commonOpts.sdp_solver_opts.verbose             = 1;

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

%% VDP

% define the problem instances based on different polynomial degrees
degV = {2  , 2:4, 1:4, 2:6};
degK = {1:2, 1:2, 1:2, 1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_VDP = generate_prob_struct_vdp_ctrl_design(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Control Design Benchmark for %s Instance %d of %d. \n', probStruct_VDP.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_vdp_ctrl_design(probStruct_VDP);

    % solve VDP example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_VDP_regularization = benchmark_fun_sq_ctrl_design(probStruct_VDP, ...
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

        solStruct_seq_VDP_hessApprox2 = benchmark_fun_sq_ctrl_design(probStruct_VDP, ...
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
    solStruct_CD_VDP = benchmark_fun_cd_ctrl_design(probStruct_VDP, ...
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
disp(cell2table(tableRows,'VariableNames',tableVariableNames));




%% Lotka

% define the problem instances based on different polynomial degrees
degV = {2  , 2:4, 1:4, 2:6};
degK = {1, 1:2, 1:2, 1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Lotka = generate_prob_struct_lotka_voltera_ctrl_design(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Control Design Benchmark for %s Instance %d of %d. \n', probStruct_Lotka.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_lotka_voltera_ctrl_design(probStruct_Lotka);

    % solve Lotka example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Lotka_regularization = benchmark_fun_sq_ctrl_design(probStruct_Lotka, ...
                                                                      commonOpts, ...
                                                                      x0, ...
                                                                      hess_est_method);

    % store in struct
    seqSolverStructs_regularization{end+1} = solStruct_seq_Lotka_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_Lotka_regularization                 = evalSequential(solStruct_seq_Lotka_regularization);
    evalSequentialStruct_Lotka_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Lotka_regularization.probNumber      = kk;

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % if enabled compute with second hessian estimation method
    if enable_second_sequential

        solStruct_seq_Lotka_hessApprox2 = benchmark_fun_sq_ctrl_design(probStruct_Lotka, ...
                                                                       commonOpts, ...
                                                                       x0, ...
                                                                       hess_est_method_2);

        % store in struct
        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_Lotka_hessApprox2;

        % just some gap in display output
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        % extract size and timing statistic from sequential
        evalSequentialStruct_Lotka_hessApprox2                 = evalSequential(solStruct_seq_Lotka_hessApprox2);
        evalSequentialStruct_Lotka_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_Lotka_hessApprox2.probNumber      = kk;
    else
        evalSequentialStruct_Lotka_hessApprox2 = [];
        solStruct_seq_Lotka_hessApprox2        = [];
    end

    % solve Lotka example with two states using custom coordinate-descent
    solStruct_CD_Lotka = benchmark_fun_cd_ctrl_design(probStruct_Lotka, ...
                                                       commonOpts, ...
                                                       x0);

    cdSolverStructs{end+1} = solStruct_CD_Lotka;
    probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_Lotka.Name, kk, numel(degV), commonOpts.solver);

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows, ...
                                probStruct_Lotka, ...
                                solStruct_CD_Lotka, ...
                                solStruct_seq_Lotka_regularization, ...
                                solStruct_seq_Lotka_hessApprox2, ...
                                true);

    % plot the scaled first order KKT condition over iterations
    if enable_plot_convergence
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_Lotka_regularization,...
                            evalSequentialStruct_Lotka_hessApprox2,...
                            solStruct_CD_Lotka)
    end
end % end of for-loop for Lotka benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames));


%% JetEngine

% define the problem instances based on different polynomial degrees
degV = {2  , 2:4, 1:4, 2:6};
degK = {1:2, 1:2, 1:2, 1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_JetEngine = generate_prob_struct_jet_engine_ctrl_design(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Control Design Benchmark for %s Instance %d of %d. \n', probStruct_JetEngine.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_jet_engine_ctrl_design(probStruct_JetEngine);

    % solve JetEngine example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_JetEngine_regularization = benchmark_fun_sq_ctrl_design(probStruct_JetEngine, ...
                                                                      commonOpts, ...
                                                                      x0, ...
                                                                      hess_est_method);

    % store in struct
    seqSolverStructs_regularization{end+1} = solStruct_seq_JetEngine_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_JetEngine_regularization                 = evalSequential(solStruct_seq_JetEngine_regularization);
    evalSequentialStruct_JetEngine_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_JetEngine_regularization.probNumber      = kk;

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % if enabled compute with second hessian estimation method
    if enable_second_sequential

        solStruct_seq_JetEngine_hessApprox2 = benchmark_fun_sq_ctrl_design(probStruct_JetEngine, ...
                                                                       commonOpts, ...
                                                                       x0, ...
                                                                       hess_est_method_2);

        % store in struct
        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_JetEngine_hessApprox2;

        % just some gap in display output
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        % extract size and timing statistic from sequential
        evalSequentialStruct_JetEngine_hessApprox2                 = evalSequential(solStruct_seq_JetEngine_hessApprox2);
        evalSequentialStruct_JetEngine_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_JetEngine_hessApprox2.probNumber      = kk;
    else
        evalSequentialStruct_JetEngine_hessApprox2 = [];
        solStruct_seq_JetEngine_hessApprox2        = [];
    end

    % solve JetEngine example with two states using custom coordinate-descent
    solStruct_CD_JetEngine = benchmark_fun_cd_ctrl_design(probStruct_JetEngine, ...
                                                           commonOpts, ...
                                                           x0);

    cdSolverStructs{end+1} = solStruct_CD_JetEngine;
    probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_JetEngine.Name, kk, numel(degV), commonOpts.solver);

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows, ...
                                probStruct_JetEngine, ...
                                solStruct_CD_JetEngine, ...
                                solStruct_seq_JetEngine_regularization, ...
                                solStruct_seq_JetEngine_hessApprox2, ...
                                true);

    % plot the scaled first order KKT condition over iterations
    if enable_plot_convergence
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_JetEngine_regularization,...
                            evalSequentialStruct_JetEngine_hessApprox2,...
                            solStruct_CD_JetEngine)
    end
end % end of for-loop for JetEngine benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))




%% MIMO

% define the problem instances based on different polynomial degrees
degV = {2  , 2:4, 2:6};
degK = {1:2, 1:2, 1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_MiMo = generate_prob_struct_mimo_ctrl_design(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Control Design Benchmark for %s Instance %d of %d. \n', probStruct_MiMo.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_mimo_ctrl_design(probStruct_MiMo);

    % solve MiMo example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_MiMo_regularization = benchmark_fun_sq_ctrl_design(probStruct_MiMo, ...
                                                                       commonOpts, ...
                                                                       x0, ...
                                                                       hess_est_method);

    % store in struct
    seqSolverStructs_regularization{end+1} = solStruct_seq_MiMo_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_MiMo_regularization                 = evalSequential(solStruct_seq_MiMo_regularization);
    evalSequentialStruct_MiMo_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_MiMo_regularization.probNumber      = kk;

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % if enabled compute with second hessian estimation method
    if enable_second_sequential

        solStruct_seq_MiMo_hessApprox2 = benchmark_fun_sq_ctrl_design(probStruct_MiMo, ...
                                                                        commonOpts, ...
                                                                        x0, ...
                                                                        hess_est_method_2);

        % store in struct
        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_MiMo_hessApprox2;

        % just some gap in display output
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        % extract size and timing statistic from sequential
        evalSequentialStruct_MiMo_hessApprox2                 = evalSequential(solStruct_seq_MiMo_hessApprox2);
        evalSequentialStruct_MiMo_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_MiMo_hessApprox2.probNumber      = kk;
    else
        evalSequentialStruct_MiMo_hessApprox2 = [];
        solStruct_seq_MiMo_hessApprox2        = [];
    end

    % solve MiMo example with two states using custom coordinate-descent
    solStruct_CD_MiMo = benchmark_fun_cd_ctrl_design(probStruct_MiMo, ...
                                                      commonOpts, ...
                                                      x0);

    cdSolverStructs{end+1} = solStruct_CD_MiMo;
    probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_MiMo.Name, kk, numel(degV), commonOpts.solver);

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows, ...
                                probStruct_MiMo, ...
                                solStruct_CD_MiMo, ...
                                solStruct_seq_MiMo_regularization, ...
                                solStruct_seq_MiMo_hessApprox2, ...
                                true);

    % plot the scaled first order KKT condition over iterations
    if enable_plot_convergence
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_MiMo_regularization,...
                            evalSequentialStruct_MiMo_hessApprox2,...
                            solStruct_CD_MiMo)
    end
end % end of for-loop for MIMO benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% DoublePendulum

% define the problem instances based on different polynomial degrees
degV = {2, 2,  2:4, 2:6};
degK = {1, 1:2,1:2, 1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_DoublePendulum = generate_prob_struct_double_pendulum_ctrl_design(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Control Design Benchmark for %s Instance %d of %d. \n', probStruct_DoublePendulum.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_double_pendulum_ctrl_design(probStruct_DoublePendulum);

    % solve DoublePendulum example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_DoublePendulum_regularization = benchmark_fun_sq_ctrl_design(probStruct_DoublePendulum, ...
                                                                       commonOpts, ...
                                                                       x0, ...
                                                                       hess_est_method);

    % store in struct
    seqSolverStructs_regularization{end+1} = solStruct_seq_DoublePendulum_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_DoublePendulum_regularization                 = evalSequential(solStruct_seq_DoublePendulum_regularization);
    evalSequentialStruct_DoublePendulum_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_DoublePendulum_regularization.probNumber      = kk;

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % if enabled compute with second hessian estimation method
    if enable_second_sequential

        solStruct_seq_DoublePendulum_hessApprox2 = benchmark_fun_sq_ctrl_design(probStruct_DoublePendulum, ...
                                                                        commonOpts, ...
                                                                        x0, ...
                                                                        hess_est_method_2);

        % store in struct
        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_DoublePendulum_hessApprox2;

        % just some gap in display output
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        % extract size and timing statistic from sequential
        evalSequentialStruct_DoublePendulum_hessApprox2                 = evalSequential(solStruct_seq_DoublePendulum_hessApprox2);
        evalSequentialStruct_DoublePendulum_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_DoublePendulum_hessApprox2.probNumber      = kk;
    else
        evalSequentialStruct_DoublePendulum_hessApprox2 = [];
        solStruct_seq_DoublePendulum_hessApprox2        = [];
    end

    % solve DoublePendulum example with two states using custom coordinate-descent
    solStruct_CD_DoublePendulum = benchmark_fun_cd_ctrl_design(probStruct_DoublePendulum, ...
                                                                commonOpts, ...
                                                                x0);

    cdSolverStructs{end+1} = solStruct_CD_DoublePendulum;
    probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_DoublePendulum.Name, kk, numel(degV), commonOpts.solver);

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows, ...
                                probStruct_DoublePendulum, ...
                                solStruct_CD_DoublePendulum, ...
                                solStruct_seq_DoublePendulum_regularization, ...
                                solStruct_seq_DoublePendulum_hessApprox2, ...
                                true);

    % plot the scaled first order KKT condition over iterations
    if enable_plot_convergence
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_DoublePendulum_regularization,...
                            evalSequentialStruct_DoublePendulum_hessApprox2,...
                            solStruct_CD_DoublePendulum)
    end
end % end of for-loop for DoublePendulum benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Aircraft (affine)

% define the problem instances based on different polynomial degrees
degV = {2  , 2:4, 1:4};
degK = {1:2, 1:2, 1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Aircraft_affine = generate_prob_struct_aircraft_affine_ctrl_design(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Control Design Benchmark for %s Instance %d of %d. \n', probStruct_Aircraft_affine.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_aircraft_affine_ctrl_design(probStruct_Aircraft_affine);

    % solve Aircraft_affine example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Aircraft_affine_regularization = benchmark_fun_sq_ctrl_design(probStruct_Aircraft_affine, ...
                                                                                  commonOpts, ...
                                                                                  x0, ...
                                                                                  hess_est_method);

    % store in struct
    seqSolverStructs_regularization{end+1} = solStruct_seq_Aircraft_affine_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_Aircraft_affine_regularization                 = evalSequential(solStruct_seq_Aircraft_affine_regularization);
    evalSequentialStruct_Aircraft_affine_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Aircraft_affine_regularization.probNumber      = kk;

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % if enabled compute with second hessian estimation method
    if enable_second_sequential

        solStruct_seq_Aircraft_affine_hessApprox2 = benchmark_fun_sq_ctrl_design(probStruct_Aircraft_affine, ...
                                                                                   commonOpts, ...
                                                                                   x0, ...
                                                                                   hess_est_method_2);

        % store in struct
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
    else
        evalSequentialStruct_Aircraft_affine_hessApprox2 = [];
        solStruct_seq_Aircraft_affine_hessApprox2        = [];
    end

    % solve Aircraft_affine example with two states using custom coordinate-descent
    solStruct_CD_Aircraft_affine = benchmark_fun_cd_ctrl_design(probStruct_Aircraft_affine, ...
                                                                 commonOpts, ...
                                                                 x0);

    cdSolverStructs{end+1} = solStruct_CD_Aircraft_affine;
    probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_Aircraft_affine.Name, kk, numel(degV), commonOpts.solver);

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows, ...
                                probStruct_Aircraft_affine, ...
                                solStruct_CD_Aircraft_affine, ...
                                solStruct_seq_Aircraft_affine_regularization, ...
                                solStruct_seq_Aircraft_affine_hessApprox2, ...
                                true);

    % plot the scaled first order KKT condition over iterations
    if enable_plot_convergence
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_Aircraft_affine_regularization,...
                            evalSequentialStruct_Aircraft_affine_hessApprox2,...
                            solStruct_CD_Aircraft_affine)
    end
end % end of for-loop for Aircraft (affine) benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Pendubot

% define the problem instances based on different polynomial degrees
degV = {2  , 2:4, 1:4};
degK = {1:2, 1:2, 1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Pendubot = generate_prob_struct_pendubot_ctrl_design(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Control Design Benchmark for %s Instance %d of %d. \n', probStruct_Pendubot.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_pendubot_ctrl_design(probStruct_Pendubot);

    % solve Pendubot example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Pendubot_regularization = benchmark_fun_sq_ctrl_design(probStruct_Pendubot, ...
                                                                           commonOpts, ...
                                                                           x0, ...
                                                                           hess_est_method);

    % store in struct
    seqSolverStructs_regularization{end+1} = solStruct_seq_Pendubot_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_Pendubot_regularization                 = evalSequential(solStruct_seq_Pendubot_regularization);
    evalSequentialStruct_Pendubot_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Pendubot_regularization.probNumber      = kk;

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % if enabled compute with second hessian estimation method
    if enable_second_sequential

        solStruct_seq_Pendubot_hessApprox2 = benchmark_fun_sq_ctrl_design(probStruct_Pendubot, ...
                                                                            commonOpts, ...
                                                                            x0, ...
                                                                            hess_est_method_2);

        % store in struct
        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_Pendubot_hessApprox2;

        % just some gap in display output
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        % extract size and timing statistic from sequential
        evalSequentialStruct_Pendubot_hessApprox2                 = evalSequential(solStruct_seq_Pendubot_hessApprox2);
        evalSequentialStruct_Pendubot_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_Pendubot_hessApprox2.probNumber      = kk;
    else
        evalSequentialStruct_Pendubot_hessApprox2 = [];
        solStruct_seq_Pendubot_hessApprox2        = [];
    end

    % solve Pendubot example with two states using custom coordinate-descent
    solStruct_CD_Pendubot = benchmark_fun_cd_ctrl_design(probStruct_Pendubot, ...
                                                          commonOpts, ...
                                                          x0);

    cdSolverStructs{end+1} = solStruct_CD_Pendubot;
    probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_Pendubot.Name, kk, numel(degV), commonOpts.solver);

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows, ...
                                probStruct_Pendubot, ...
                                solStruct_CD_Pendubot, ...
                                solStruct_seq_Pendubot_regularization, ...
                                solStruct_seq_Pendubot_hessApprox2, ...
                                true);

    % plot the scaled first order KKT condition over iterations
    if enable_plot_convergence
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_Pendubot_regularization,...
                            evalSequentialStruct_Pendubot_hessApprox2,...
                            solStruct_CD_Pendubot)
    end
end % end of for-loop for Pendubot benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Spacecraft_orbit

% define the problem instances based on different polynomial degrees
degV = {2,   2, 2:4};
degK = {1, 1:2, 1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Spacecraft_orbit = generate_prob_struct_spacecraft_orbit_ctrl_design(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Control Design Benchmark for %s Instance %d of %d. \n', probStruct_Spacecraft_orbit.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_spacecraft_orbit_ctrl_design(probStruct_Spacecraft_orbit);

    % solve Spacecraft_orbit example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Spacecraft_orbit_regularization = benchmark_fun_sq_ctrl_design(probStruct_Spacecraft_orbit, ...
                                                                            commonOpts, ...
                                                                            x0, ...
                                                                            hess_est_method);

    % store in struct
    seqSolverStructs_regularization{end+1} = solStruct_seq_Spacecraft_orbit_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_Spacecraft_orbit_regularization                 = evalSequential(solStruct_seq_Spacecraft_orbit_regularization);
    evalSequentialStruct_Spacecraft_orbit_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Spacecraft_orbit_regularization.probNumber      = kk;

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % if enabled compute with second hessian estimation method
    if enable_second_sequential

        solStruct_seq_Spacecraft_orbit_hessApprox2 = benchmark_fun_sq_ctrl_design(probStruct_Spacecraft_orbit, ...
                                                                             commonOpts, ...
                                                                             x0, ...
                                                                             hess_est_method_2);

        % store in struct
        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_Spacecraft_orbit_hessApprox2;

        % just some gap in display output
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        % extract size and timing statistic from sequential
        evalSequentialStruct_Spacecraft_orbit_hessApprox2                 = evalSequential(solStruct_seq_Spacecraft_orbit_hessApprox2);
        evalSequentialStruct_Spacecraft_orbit_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_Spacecraft_orbit_hessApprox2.probNumber      = kk;
    else
        evalSequentialStruct_Spacecraft_orbit_hessApprox2 = [];
        solStruct_seq_Spacecraft_orbit_hessApprox2        = [];
    end

    % solve Spacecraft_orbit example with two states using custom coordinate-descent
    solStruct_CD_Spacecraft_orbit = benchmark_fun_cd_ctrl_design(probStruct_Spacecraft_orbit, ...
                                                                  commonOpts, ...
                                                                  x0);

    cdSolverStructs{end+1} = solStruct_CD_Spacecraft_orbit;
    probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_Spacecraft_orbit.Name, kk, numel(degV), commonOpts.solver);

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows, ...
                                probStruct_Spacecraft_orbit, ...
                                solStruct_CD_Spacecraft_orbit, ...
                                solStruct_seq_Spacecraft_orbit_regularization, ...
                                solStruct_seq_Spacecraft_orbit_hessApprox2, ...
                                true);

    % plot the scaled first order KKT condition over iterations
    if enable_plot_convergence
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_Spacecraft_orbit_regularization,...
                            evalSequentialStruct_Spacecraft_orbit_hessApprox2,...
                            solStruct_CD_Spacecraft_orbit)
    end
end % end of for-loop for Spacecraft_orbit benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Satellite

% define the problem instances based on different polynomial degrees
degV = {2,   2, 2:4};
degK = {1, 1:2, 1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Satellite = generate_prob_struct_satellite_ctrl_design(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Control Design Benchmark for %s Instance %d of %d. \n', probStruct_Satellite.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_satellite_ctrl_design(probStruct_Satellite);

    % solve Satellite example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Satellite_regularization = benchmark_fun_sq_ctrl_design(probStruct_Satellite, ...
                                                                            commonOpts, ...
                                                                            x0, ...
                                                                            hess_est_method);

    % store in struct
    seqSolverStructs_regularization{end+1} = solStruct_seq_Satellite_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_Satellite_regularization                 = evalSequential(solStruct_seq_Satellite_regularization);
    evalSequentialStruct_Satellite_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Satellite_regularization.probNumber      = kk;

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % if enabled compute with second hessian estimation method
    if enable_second_sequential

        solStruct_seq_Satellite_hessApprox2 = benchmark_fun_sq_ctrl_design(probStruct_Satellite, ...
                                                                             commonOpts, ...
                                                                             x0, ...
                                                                             hess_est_method_2);

        % store in struct
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
    else
        evalSequentialStruct_Satellite_hessApprox2 = [];
        solStruct_seq_Satellite_hessApprox2        = [];
    end

    % solve Satellite example with two states using custom coordinate-descent
    solStruct_CD_Satellite = benchmark_fun_cd_ctrl_design(probStruct_Satellite, ...
                                                           commonOpts, ...
                                                           x0);

    cdSolverStructs{end+1} = solStruct_CD_Satellite;
    probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_Satellite.Name, kk, numel(degV), commonOpts.solver);

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
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

%% Evaluation


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

% evaluate benchmark
plotPerformanceProfiles(seqMethodStructs, ...
                        cdSolverStructs, ...
                        solverNames, ...
                        commonOpts.maxSolveTime, ...
                        probNames)

% final results table accumulated over all runs and all problems
resultsTable = cell2table(tableRows,'VariableNames',tableVariableNames);
disp(resultsTable)

% remove the added folders from the path again
for ii = 1:numel(addPaths)
    rmpath(addPaths{ii})
end
