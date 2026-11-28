%% ***********************************************************************
%
%  Region-of-attraction (ROA) Estimation Benchmark Runner
%
%
%   SOS Problem for ROA Estimation
%
%   min      ||g - (V-β)||
%   V,s 
%
%   s.j.    s                        is SOS
%           V  - ϵ(x)                is SOS
%           s(V-β) - ∇Vf(x,u) - ϵ(x) is SOS
%
%   V    - Lyapunov function
%   s    - SOS multiplier
%   β    - Level set (fixed)
%   ϵ(x) - Small term
%   g    - Safe set/desired operating region
%  
%
%  Overview
%  
%  Benchmark 1: VDP
%  Benchmark 2: Multi-Input-Multi-Output
%  Benchmark 3: VDP_coupled
%  Benchmark 4: Double Pendulum
%  Benchmark 5: Aircraft 4D
%  Benchmark 6: Satellite
%  Benchmark 7: F/A-18
%  Benchmark 8: 12-State Dynamical System
%
%
%  License: GNU GENERAL PUBLIC LICENSE Version 3
%
% ***********************************************************************
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
name_of_benchmark = 'ROA_benchmark';

% get versions
matlabVersion = version;
casadiVersion = casadi.CasadiMeta.version;


% add the problem specific folders and the shared helper functions
runnerDir = pwd;
addPaths  = { ...
    fullfile(runnerDir,'vdp_roa'), ...
    fullfile(runnerDir,'vdp_coupled_roa'), ...
    fullfile(runnerDir,'mimo_roa'), ...
    fullfile(runnerDir,'double_pendulum_roa'), ...
    fullfile(runnerDir,'aircraft_roa'), ...
    fullfile(runnerDir,'spacecraft_att_roa'), ...
    fullfile(runnerDir,'f18_roa'), ...
    fullfile(runnerDir,'network_roa'), ...
    fullfile(runnerDir,'..','helperFunctions')};

for ii = 1:numel(addPaths)
    addpath(addPaths{ii})
end

% enable plotting of convergence for each run
enable_plot_convergence  = true;

% enable runs with second hessian estimation method
enable_second_sequential = true;

% Hessian estimation methods of the sequential solver to compare
hess_est_method_1 = 'regularization_frobenius'; % our default
hess_est_method_2 = 'regularization_levenberg'; % 'damped_bfgs', 'regularization_levenberg', 'regularization_mirror', 'regularization_frobenius'

% common options shared by the methods
commonOpts.maxiter          = 500;
commonOpts.tolerance_opt    = 1e-4;
commonOpts.eps_cost         = 1e-6;
commonOpts.small_prog_count = 10;
commonOpts.maxSolveTime     = 3600;
commonOpts.almost_opt_count = 10;
commonOpts.solver           = 'mosek'; % 'clarabel', 'scs'
commonOpts.miniter          = 5; 

% SDP solver settings
if strcmpi(commonOpts.solver,'mosek')

    commonOpts.sdp_solver_opts.mosek_param = struct('MSK_IPAR_NUM_THREADS',1);
    commonOpts.sdp_solver_opts.mosek_echo = 0;

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
% (if selected)
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
probNames                       = {};


total_benchmark_time_start = tic;

%% VDP

% define the problem instances based on different polynomial degrees
degV = {2,2:4,1:4,2:6};
degS = {2,2,2:4,2:4};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_VDP = generate_prob_struct_vdp_roa(degV{kk},degS{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_VDP.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_vdp_roa(probStruct_VDP);

    % solve VDP example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_VDP_regularization = benchmark_fun_sq_roa(probStruct_VDP, ...
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

        solStruct_seq_VDP_hessApprox2 = benchmark_fun_sq_roa(probStruct_VDP, ...
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

    % solve VDP example using custom coordinate-descent
    solStruct_CD_VDP = benchmark_fun_cd_roa(probStruct_VDP, ...
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



%% Multi-Input-Multi-Output

% define the problem instances based on different polynomial degrees
degV = {2, 2:4, 2:4, 2:6};
degS = {2,   2, 2:4, 2:4};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_MiMo = generate_prob_struct_mimo_roa(degV{kk},degS{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_MiMo.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_mimo_roa(probStruct_MiMo);

    % solve MiMo example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_MiMo_regularization = benchmark_fun_sq_roa(probStruct_MiMo, ...
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

        solStruct_seq_MiMo_hessApprox2 = benchmark_fun_sq_roa(probStruct_MiMo, ...
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
    solStruct_CD_MiMo = benchmark_fun_cd_roa(probStruct_MiMo, ...
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

end % end of for-loop for Multi-Input-Multi-Output benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))




%% VDP_coupled

% define the problem instances based on different polynomial degrees
degV = {2,2:4,1:4,2:6};
degS = {2,2,2:4,2:4};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_VDP_coupled = generate_prob_struct_vdp_coupled_roa(degV{kk},degS{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_VDP_coupled.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_vdp_coupled_roa(probStruct_VDP_coupled);

    % solve VDP_coupled example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_VDP_coupled_regularization = benchmark_fun_sq_roa(probStruct_VDP_coupled, ...
                                                               commonOpts, ...
                                                               x0, ...
                                                               hess_est_method);

    % store in struct
    seqSolverStructs_regularization{end+1} = solStruct_seq_VDP_coupled_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_VDP_coupled_regularization                 = evalSequential(solStruct_seq_VDP_coupled_regularization);
    evalSequentialStruct_VDP_coupled_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_VDP_coupled_regularization.probNumber      = kk;

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % if enabled compute with second hessian estimation method
    if enable_second_sequential

        solStruct_seq_VDP_coupled_hessApprox2 = benchmark_fun_sq_roa(probStruct_VDP_coupled, ...
                                                                commonOpts, ...
                                                                x0, ...
                                                                hess_est_method_2);

        % store in struct
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
    else
        evalSequentialStruct_VDP_coupled_hessApprox2 = [];
        solStruct_seq_VDP_coupled_hessApprox2        = [];
    end

    % solve VDP_coupled example with two states using custom coordinate-descent
    solStruct_CD_VDP_coupled = benchmark_fun_cd_roa(probStruct_VDP_coupled, ...
                                                     commonOpts, ...
                                                     x0);

    cdSolverStructs{end+1} = solStruct_CD_VDP_coupled;
    probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_VDP_coupled.Name, kk, numel(degV), commonOpts.solver);

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows, ...
                                probStruct_VDP_coupled, ...
                                solStruct_CD_VDP_coupled, ...
                                solStruct_seq_VDP_coupled_regularization, ...
                                solStruct_seq_VDP_coupled_hessApprox2, ...
                                true);

    % plot the scaled first order KKT condition over iterations
    if enable_plot_convergence
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_VDP_coupled_regularization,...
                            evalSequentialStruct_VDP_coupled_hessApprox2,...
                            solStruct_CD_VDP_coupled)
    end

end % end of for-loop for VDP_coupled benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Double Pendulum

% define the problem instances based on different polynomial degrees
degV = {2;2:4;2:4};
degS = {2;2;2:4};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_DoublePendulum = generate_prob_struct_double_pendulum_roa(degV{kk},degS{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_DoublePendulum.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_double_pendulum_roa(probStruct_DoublePendulum);

    % solve DoublePendulum example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_DoublePendulum_regularization = benchmark_fun_sq_roa(probStruct_DoublePendulum, ...
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

        solStruct_seq_DoublePendulum_hessApprox2 = benchmark_fun_sq_roa(probStruct_DoublePendulum, ...
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
    solStruct_CD_DoublePendulum = benchmark_fun_cd_roa(probStruct_DoublePendulum, ...
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

end % end of for-loop for Double Pendulum benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Aircraft 4D

% define the problem instances based on different polynomial degrees
degV = {2,   2, 2:4, 1:4};
degS = {2, 2:4, 2:4, 2:4};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Aircraft = generate_prob_struct_aircraft_roa(degV{kk},degS{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_Aircraft.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_aircraft_roa(probStruct_Aircraft);

    % solve Aircraft example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Aircraft_regularization = benchmark_fun_sq_roa(probStruct_Aircraft, ...
                                                                    commonOpts, ...
                                                                    x0, ...
                                                                    hess_est_method);

    % store in struct
    seqSolverStructs_regularization{end+1} = solStruct_seq_Aircraft_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_Aircraft_regularization                 = evalSequential(solStruct_seq_Aircraft_regularization);
    evalSequentialStruct_Aircraft_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Aircraft_regularization.probNumber      = kk;

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % if enabled compute with second hessian estimation method
    if enable_second_sequential

        solStruct_seq_Aircraft_hessApprox2 = benchmark_fun_sq_roa(probStruct_Aircraft, ...
                                                                     commonOpts, ...
                                                                     x0, ...
                                                                     hess_est_method_2);

        % store in struct
        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_Aircraft_hessApprox2;

        % just some gap in display output
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        % extract size and timing statistic from sequential
        evalSequentialStruct_Aircraft_hessApprox2                 = evalSequential(solStruct_seq_Aircraft_hessApprox2);
        evalSequentialStruct_Aircraft_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_Aircraft_hessApprox2.probNumber      = kk;
    else
        evalSequentialStruct_Aircraft_hessApprox2 = [];
        solStruct_seq_Aircraft_hessApprox2        = [];
    end

    % solve Aircraft example with two states using custom coordinate-descent
    solStruct_CD_Aircraft = benchmark_fun_cd_roa(probStruct_Aircraft, ...
                                                  commonOpts, ...
                                                  x0);

    cdSolverStructs{end+1} = solStruct_CD_Aircraft;
    probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_Aircraft.Name, kk, numel(degV), commonOpts.solver);

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
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

end % end of for-loop for Aircraft 4D benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))



%% Satellite

% define the problem instances based on different polynomial degrees
degV = {2, 2:4, 2:4, 1:4};
degS = {2,   2, 2:4, 2:4};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Satellite = generate_prob_struct_spacecraft_att_roa(degV{kk},degS{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_Satellite.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_spacecraft_att_roa(probStruct_Satellite);

    % solve Satellite example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Satellite_regularization = benchmark_fun_sq_roa(probStruct_Satellite, ...
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

        solStruct_seq_Satellite_hessApprox2 = benchmark_fun_sq_roa(probStruct_Satellite, ...
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
    solStruct_CD_Satellite = benchmark_fun_cd_roa(probStruct_Satellite, ...
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



%% F/A-18

% define the problem instances based on different polynomial degrees
degV = {2, 2:4, 1:4};
degS = {2,   2, 2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_F18 = generate_prob_struct_f18_roa(degV{kk},degS{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_F18.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_init_guess_f18_roa(probStruct_F18);

    % solve F-18 example with seven states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_7D_regularization = benchmark_fun_sq_roa(probStruct_F18, ...
                                                              commonOpts, ...
                                                              x0, ...
                                                              hess_est_method);

    % store in struct
    seqSolverStructs_regularization{end+1} = solStruct_seq_7D_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_7D_regularization                 = evalSequential(solStruct_seq_7D_regularization);
    evalSequentialStruct_7D_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_7D_regularization.probNumber      = kk;

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % if enabled compute with second hessian estimation method
    if enable_second_sequential

        solStruct_seq_7D_hessApprox2 = benchmark_fun_sq_roa(probStruct_F18, ...
                                                               commonOpts, ...
                                                               x0, ...
                                                               hess_est_method_2);

        % store in struct
        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_7D_hessApprox2;

        % just some gap in display output
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')
        fprintf('\n')

        % extract size and timing statistic from sequential
        evalSequentialStruct_7D_hessApprox2                 = evalSequential(solStruct_seq_7D_hessApprox2);
        evalSequentialStruct_7D_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_7D_hessApprox2.probNumber      = kk;
    else
        evalSequentialStruct_7D_hessApprox2 = [];
        solStruct_seq_7D_hessApprox2        = [];
    end

    % solve F-18 example with seven states using custom coordinate-descent
    solStruct_CD_7D = benchmark_fun_cd_roa(probStruct_F18, ...
                                             commonOpts, ...
                                             x0);

    cdSolverStructs{end+1} = solStruct_CD_7D;
    probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_F18.Name, kk, numel(degV), commonOpts.solver);

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows, ...
                                probStruct_F18, ...
                                solStruct_CD_7D, ...
                                solStruct_seq_7D_regularization, ...
                                solStruct_seq_7D_hessApprox2, ...
                                true);

    % plot the scaled first order KKT condition over iterations
    if enable_plot_convergence
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_7D_regularization,...
                            evalSequentialStruct_7D_hessApprox2,...
                            solStruct_CD_7D)
    end

end % end of for-loop for F/A-18 benchmark

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))

%% 12-State Dynamical System

% generate problem structure (dynamics, polynomial dec. variables)
probStruct_12State = generate_prob_struct_12state_roa();
fprintf('------------------------------------------------------------------------------------------ \n')
fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_12State.Name, 1, 1)
fprintf('------------------------------------------------------------------------------------------\n')

% generate an initial guess
x0 = generate_init_guess_12state_roa(probStruct_12State);

% solve 12 state dynamical system example with 12 states using sequential SOS
hess_est_method = hess_est_method_1;
solStruct_seq_12D_regularization = benchmark_fun_sq_roa(probStruct_12State, ...
                                                           commonOpts, ...
                                                           x0, ...
                                                           hess_est_method);

% store in struct
seqSolverStructs_regularization{end+1} = solStruct_seq_12D_regularization;

% extract size and timing statistic from sequential and add some meta
% data for plotting
evalSequentialStruct_12D_regularization                 = evalSequential(solStruct_seq_12D_regularization);
evalSequentialStruct_12D_regularization.hess_est_method = hess_est_method_1;
evalSequentialStruct_12D_regularization.probNumber      = 1;

% just some gap in display output
fprintf('\n')
fprintf('\n')
fprintf('\n')
fprintf('\n')

% if enabled compute with second hessian estimation method
if enable_second_sequential

    solStruct_seq_12D_hessApprox2 = benchmark_fun_sq_roa(probStruct_12State, ...
                                                            commonOpts, ...
                                                            x0, ...
                                                            hess_est_method_2);

    % store in struct
    seqSolverStructs_hessApprox2{end+1} = solStruct_seq_12D_hessApprox2;

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % extract size and timing statistic from sequential
    evalSequentialStruct_12D_hessApprox2                 = evalSequential(solStruct_seq_12D_hessApprox2);
    evalSequentialStruct_12D_hessApprox2.hess_est_method = hess_est_method_2;
    evalSequentialStruct_12D_hessApprox2.probNumber      = 1;
else
    evalSequentialStruct_12D_hessApprox2 = [];
    solStruct_seq_12D_hessApprox2        = [];
end

% solve 12 state dynamical system example with 12 states using custom coordinate-descent
solStruct_CD_12D = benchmark_fun_cd_roa(probStruct_12State, ...
                                         commonOpts, ...
                                         x0);

cdSolverStructs{end+1} = solStruct_CD_12D;
probNames{end+1}       = sprintf('%s (instance %d of %d with %s)', probStruct_12State.Name, 1, 1, commonOpts.solver);

% just some gap in display output
fprintf('\n')
fprintf('\n')
fprintf('\n')
fprintf('\n')

% append the current run to the accumulating results table
tableRows = appendResultRow(tableRows, ...
                            probStruct_12State, ...
                            solStruct_CD_12D, ...
                            solStruct_seq_12D_regularization, ...
                            solStruct_seq_12D_hessApprox2, ...
                            true);

% plot the scaled first order KKT condition over iterations
if enable_plot_convergence
        plotConvergence(probNames{end}, ...
                        enable_second_sequential,...
                        evalSequentialStruct_12D_regularization,...
                        evalSequentialStruct_12D_hessApprox2,...
                        solStruct_CD_12D)
end

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
% are done (exactly the folders added at the top of this script)
for ii = 1:numel(addPaths)
    rmpath(addPaths{ii})
end
