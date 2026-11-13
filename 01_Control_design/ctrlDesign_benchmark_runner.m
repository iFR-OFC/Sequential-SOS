%% ***********************************************************************
%
%  Control Design Benchmark Test Runner
%
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
    fullfile(runnerDir,'VDP'), ...
    fullfile(runnerDir,'MultiInput'), ...
    fullfile(runnerDir,'Aircraft_affine'), ...
    fullfile(runnerDir,'Pendubot'),...
    fullfile(runnerDir,'Quadcopter'),...
    fullfile(runnerDir,'Satellite'), ...
    fullfile(runnerDir,'..','helperFunctions')};

for ii = 1:numel(addPaths)
    addpath(addPaths{ii})
end

% enable plotting of convergence for each run
enable_plot_convergence  = true; 

% enable second run of sequential SOS with different Hessian estimation
enable_second_sequential = true;

% Hessian estimation methods of the sequential solver to compare
hess_est_method_1 = 'regularization'; % our default (see also sequentialCommon.m)
hess_est_method_2 =  'BFGS';          % 'BFGS', 'levenberg', 'mirroring'

% common options shared by the methods
commonOpts.maxiter       = 500;
commonOpts.tolerance_opt = 1e-4;
commonOpts.maxSolveTime  = 3600;
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

    commonOpts.sdp_solver_opts.dense               = 1;
    commonOpts.sdp_solver_opts.verbose             = 0;

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


%% VDP

% define the problem instances based on different polynomial degrees
degV = {2  , 2:4, 1:4, 2:6};
degK = {1:2, 1:2, 1:2, 1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_VDP = generateVDP_probStruct(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Control Design Benchmark for %s Instance %d of %d. \n', probStruct_VDP.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_VDP(probStruct_VDP);

    % solve VDP example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_VDP_regularization = solve_ctrlDesign_seq_bench_fun(probStruct_VDP, ...
                                                                      commonOpts, ...
                                                                      x0, ...
                                                                      hess_est_method);

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

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_VDP_hessApprox2 = solve_ctrlDesign_seq_bench_fun(probStruct_VDP, ...
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

    % solve VDP example with two states using custom coordinate-descent
    solStruct_CD_VDP = solve_ctrlDesign_CD_bench_fun(probStruct_VDP,commonOpts,x0);

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


%% MIMO

% define the problem instances based on different polynomial degrees
degV = {2  , 2:4, 2:6};
degK = {1:2, 1:2, 1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_MiMo = generateMiMo_probStruct(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Control Design Benchmark for %s Instance %d of %d. \n', probStruct_MiMo.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_MiMo(probStruct_MiMo);

    % solve MiMo example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_MiMo_regularization = solve_ctrlDesign_seq_bench_fun(probStruct_MiMo, ...
                                                                       commonOpts, ...
                                                                       x0, ...
                                                                       hess_est_method);

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

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_MiMo_hessApprox2 = solve_ctrlDesign_seq_bench_fun(probStruct_MiMo, ...
                                                                        commonOpts, ...
                                                                        x0, ...
                                                                        hess_est_method);

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
    end

    % solve MiMo example with two states using custom coordinate-descent
    solStruct_CD_MiMo = solve_ctrlDesign_CD_bench_fun(probStruct_MiMo,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_MiMo;
    probNames{end+1}       = sprintf('%s (instance %d of %d)', probStruct_MiMo.Name, kk, numel(degV));

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows, ...
                                    probStruct_MiMo, ...
                                    solStruct_CD_MiMo, ...
                                    solStruct_seq_MiMo_regularization, ...
                                    solStruct_seq_MiMo_hessApprox2, ...
                                    true);
    else
        tableRows = appendResultRow(tableRows, ...
                                    probStruct_MiMo, ...
                                    solStruct_CD_MiMo, ...
                                    solStruct_seq_MiMo_regularization, ...
                                    [], ...
                                    true);
    end

    if enable_plot_convergence
        if enable_second_sequential
            % Evaluate
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_MiMo_regularization,...
                            evalSequentialStruct_MiMo_hessApprox2,...
                            solStruct_CD_MiMo)
        else
            % Evaluate
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_MiMo_regularization,...
                            [],...
                            solStruct_CD_MiMo)
        end
    end
end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Aircraft (affine)

% define the problem instances based on different polynomial degrees
degV = {2  , 2:4, 1:4};
degK = {1:2, 1:2, 1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Aircraft_affine = generateAircraft_affine_probStruct(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Control Design Benchmark for %s Instance %d of %d. \n', probStruct_Aircraft_affine.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_Aircraft_affine(probStruct_Aircraft_affine);

    % solve Aircraft_affine example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Aircraft_affine_regularization = solve_ctrlDesign_seq_bench_fun(probStruct_Aircraft_affine, ...
                                                                                  commonOpts, ...
                                                                                  x0, ...
                                                                                  hess_est_method);

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

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_Aircraft_affine_hessApprox2 = solve_ctrlDesign_seq_bench_fun(probStruct_Aircraft_affine, ...
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

    % solve Aircraft_affine example with two states using custom coordinate-descent
    solStruct_CD_Aircraft_affine = solve_ctrlDesign_CD_bench_fun(probStruct_Aircraft_affine,commonOpts,x0);

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


%% Pendubot

% define the problem instances based on different polynomial degrees
degV = {2  , 2:4, 1:4};
degK = {1:2, 1:2, 1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Pendubot = generatePendubot_probStruct(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Control Design Benchmark for %s Instance %d of %d. \n', probStruct_Pendubot.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_Pendubot(probStruct_Pendubot);

    % solve Pendubot example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Pendubot_regularization = solve_ctrlDesign_seq_bench_fun(probStruct_Pendubot, ...
                                                                           commonOpts, ...
                                                                           x0, ...
                                                                           hess_est_method);

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

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_Pendubot_hessApprox2 = solve_ctrlDesign_seq_bench_fun(probStruct_Pendubot, ...
                                                                            commonOpts, ...
                                                                            x0, ...
                                                                            hess_est_method);

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
    end

    % solve Pendubot example with two states using custom coordinate-descent
    solStruct_CD_Pendubot = solve_ctrlDesign_CD_bench_fun(probStruct_Pendubot,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_Pendubot;
    probNames{end+1}       = sprintf('%s (instance %d of %d)', probStruct_Pendubot.Name, kk, numel(degV));

    % just some gap in display output
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')
    fprintf('\n')

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows, ...
                                    probStruct_Pendubot, ...
                                    solStruct_CD_Pendubot, ...
                                    solStruct_seq_Pendubot_regularization, ...
                                    solStruct_seq_Pendubot_hessApprox2, ...
                                    true);
    else
        tableRows = appendResultRow(tableRows, ...
                                    probStruct_Pendubot, ...
                                    solStruct_CD_Pendubot, ...
                                    solStruct_seq_Pendubot_regularization, ...
                                    [], ...
                                    true);
    end

    if enable_plot_convergence
        if enable_second_sequential
            % Evaluate
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_Pendubot_regularization,...
                            evalSequentialStruct_Pendubot_hessApprox2,...
                            solStruct_CD_Pendubot)
        else
            % Evaluate
            plotConvergence(probNames{end}, ...
                            enable_second_sequential,...
                            evalSequentialStruct_Pendubot_regularization,...
                            [],...
                            solStruct_CD_Pendubot)
        end
    end
end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))

%% Satellite

% define the problem instances based on different polynomial degrees
degV = {2,   2, 2:4};
degK = {1, 1:2, 1:2};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Satellite = generateSatellite_probStruct(degV{kk},degK{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Control Design Benchmark for %s Instance %d of %d. \n', probStruct_Satellite.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_Satellite(probStruct_Satellite);

    % solve Satellite example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Satellite_regularization = solve_ctrlDesign_seq_bench_fun(probStruct_Satellite, ...
                                                                            commonOpts, ...
                                                                            x0, ...
                                                                            hess_est_method);

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

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_Satellite_hessApprox2 = solve_ctrlDesign_seq_bench_fun(probStruct_Satellite, ...
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

    % solve Satellite example with two states using custom coordinate-descent
    solStruct_CD_Satellite = solve_ctrlDesign_CD_bench_fun(probStruct_Satellite,commonOpts,x0);

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


%% Store benchmark data

fprintf('\n')
fprintf('Store benchmark data in .mat file\n')
fprintf('\n')

save([name_of_benchmark '_' commonOpts.solver '.mat'],...         % name of .mat files
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
