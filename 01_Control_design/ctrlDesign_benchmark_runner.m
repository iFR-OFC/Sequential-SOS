%% ***********************************************************************
%
%  Control Design Benchmark Test Runner
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

% add the problem specific folders
addpath("VDP/")
addpath("MultiInput/")
addpath("Aircraft_affine/")
addpath("DoublePendulum/")
addpath("Satellite/")

% add the shared helper functions (evaluation and plotting)
addpath("../helperFunctions/")

% enable plotting of convergence for each run
enable_plot_convergence = false;
enable_second_sequential = false;

% Hessian estimation methods of the sequential solver to compare
hess_est_method_1 = 'regularization'; % our default (see also sequentialCommon.m)
hess_est_method_2 = 'mirroring';      % 'BFGS', 'levenberg', 'mirroring'

% common options shared by the methods
commonOpts.maxiter       = 350;
commonOpts.tolerance_opt = 1e-4;
commonOpts.maxSolveTime  = 3600;
commonOpts.solver        = 'mosek';


% variable names of the results table (one row per run)
tableVariableNames = {'Problem Name','no. con. (SDP)','no. dec. (SDP)','Size Gram matrices', ...
    'Build Time (CD)','Build Time (SQ)','Solve Time (CD)','Solve Time (SQ)', ...
    'Iterations (CD)','Iterations (SQ)','Final Cost (CD)','Final Cost (SQ)'};

% extend the results table by the sequential run with the second Hessian
% estimation method (cf. enable_second_sequential)
if enable_second_sequential
    tableVariableNames = [tableVariableNames, ...
        ['Build Time (SQ ' hess_est_method_2 ')'], ...
        ['Solve Time (SQ ' hess_est_method_2 ')'], ...
        ['Iterations (SQ ' hess_est_method_2 ')'], ...
        ['Final Cost (SQ ' hess_est_method_2 ')']];
end

% initialize the cell array that accumulates one row per run of the
% results table (one column per table variable)
tableRows = cell(0, numel(tableVariableNames));

% initialize the cell arrays of stored runs (one entry per run)
seqSolverStructs_regularization = {};
seqSolverStructs_hessApprox2    = {};
cdSolverStructs                 = {};


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
    solStruct_seq_VDP_regularization = solve_ctrlDesign_seq_bench_fun(probStruct_VDP,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_VDP_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_VDP_regularization                 = evalSequential(solStruct_seq_VDP_regularization);
    evalSequentialStruct_VDP_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_VDP_regularization.probNumber      = kk;

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_VDP_hessApprox2 = solve_ctrlDesign_seq_bench_fun(probStruct_VDP,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_VDP_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_VDP_hessApprox2 = evalSequential(solStruct_seq_VDP_hessApprox2);
        evalSequentialStruct_VDP_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_VDP_hessApprox2.probNumber = kk;
    end

    % solve VDP example with two states using custom coordinate-descent
    solStruct_CD_VDP = solve_ctrlDesign_CD_bench_fun(probStruct_VDP,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_VDP;

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows,probStruct_VDP,solStruct_CD_VDP, ...
            solStruct_seq_VDP_regularization,solStruct_seq_VDP_hessApprox2);
    else
        tableRows = appendResultRow(tableRows,probStruct_VDP,solStruct_CD_VDP, ...
                                    solStruct_seq_VDP_regularization);
    end

    if enable_plot_convergence
        % Evaluate
        plotConvergence(probStruct_VDP, ...
                        enable_second_sequential,...
                        evalSequentialStruct_VDP_regularization,...
                        evalSequentialStruct_VDP_hessApprox2,...
                        solStruct_CD_VDP)
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
    solStruct_seq_MiMo_regularization = solve_ctrlDesign_seq_bench_fun(probStruct_MiMo,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_MiMo_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_MiMo_regularization                 = evalSequential(solStruct_seq_MiMo_regularization);
    evalSequentialStruct_MiMo_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_MiMo_regularization.probNumber      = kk;

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_MiMo_hessApprox2 = solve_ctrlDesign_seq_bench_fun(probStruct_MiMo,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_MiMo_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_MiMo_hessApprox2 = evalSequential(solStruct_seq_MiMo_hessApprox2);
        evalSequentialStruct_MiMo_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_MiMo_hessApprox2.probNumber = kk;
    end

    % solve MiMo example with two states using custom coordinate-descent
    solStruct_CD_MiMo = solve_ctrlDesign_CD_bench_fun(probStruct_MiMo,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_MiMo;

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows,probStruct_MiMo,solStruct_CD_MiMo, ...
            solStruct_seq_MiMo_regularization,solStruct_seq_MiMo_hessApprox2);
    else
        tableRows = appendResultRow(tableRows,probStruct_MiMo,solStruct_CD_MiMo, ...
            solStruct_seq_MiMo_regularization);
    end

    if enable_plot_convergence
        % Evaluate
        plotConvergence(probStruct_MiMo, ...
            enable_second_sequential,...
            evalSequentialStruct_MiMo_regularization,...
            evalSequentialStruct_MiMo_hessApprox2,...
            solStruct_CD_MiMo)
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
    solStruct_seq_Aircraft_affine_regularization = solve_ctrlDesign_seq_bench_fun(probStruct_Aircraft_affine,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_Aircraft_affine_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_Aircraft_affine_regularization                 = evalSequential(solStruct_seq_Aircraft_affine_regularization);
    evalSequentialStruct_Aircraft_affine_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Aircraft_affine_regularization.probNumber      = kk;

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_Aircraft_affine_hessApprox2 = solve_ctrlDesign_seq_bench_fun(probStruct_Aircraft_affine,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_Aircraft_affine_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_Aircraft_affine_hessApprox2 = evalSequential(solStruct_seq_Aircraft_affine_hessApprox2);
        evalSequentialStruct_Aircraft_affine_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_Aircraft_affine_hessApprox2.probNumber = kk;
    end

    % solve Aircraft_affine example with two states using custom coordinate-descent
    solStruct_CD_Aircraft_affine = solve_ctrlDesign_CD_bench_fun(probStruct_Aircraft_affine,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_Aircraft_affine;

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows,probStruct_Aircraft_affine,solStruct_CD_Aircraft_affine, ...
            solStruct_seq_Aircraft_affine_regularization,solStruct_seq_Aircraft_affine_hessApprox2);
    else
        tableRows = appendResultRow(tableRows,probStruct_Aircraft_affine,solStruct_CD_Aircraft_affine, ...
            solStruct_seq_Aircraft_affine_regularization);
    end

    if enable_plot_convergence
        % Evaluate
        plotConvergence(probStruct_Aircraft_affine, ...
            enable_second_sequential,...
            evalSequentialStruct_Aircraft_affine_regularization,...
            evalSequentialStruct_Aircraft_affine_hessApprox2,...
            solStruct_CD_Aircraft_affine)
    end
end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))

%% Double Pendulum

% define the problem instances based on different polynomial degrees
degV = {2  , 2:4, 2:4};
degK = {1, 1, 1:2};

for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_DoublePendulum = generateDoublePendulum_probStruct(degV{kk},degK{kk});
    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Control Design Benchmark for %s Instance %d of %d. \n', probStruct_DoublePendulum.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_DoublePendulum(probStruct_DoublePendulum);

    % solve DoublePendulum example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_DoublePendulum_regularization = solve_ctrlDesign_seq_bench_fun(probStruct_DoublePendulum,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_DoublePendulum_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_DoublePendulum_regularization                 = evalSequential(solStruct_seq_DoublePendulum_regularization);
    evalSequentialStruct_DoublePendulum_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_DoublePendulum_regularization.probNumber      = kk;

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_DoublePendulum_hessApprox2 = solve_ctrlDesign_seq_bench_fun(probStruct_DoublePendulum,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_DoublePendulum_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_DoublePendulum_hessApprox2 = evalSequential(solStruct_seq_DoublePendulum_hessApprox2);
        evalSequentialStruct_DoublePendulum_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_DoublePendulum_hessApprox2.probNumber = kk;
    end

    % solve DoublePendulum example with two states using custom coordinate-descent
    solStruct_CD_DoublePendulum = solve_ctrlDesign_CD_bench_fun(probStruct_DoublePendulum,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_DoublePendulum;

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows,probStruct_DoublePendulum,solStruct_CD_DoublePendulum, ...
            solStruct_seq_DoublePendulum_regularization,solStruct_seq_DoublePendulum_hessApprox2);
    else
        tableRows = appendResultRow(tableRows,probStruct_DoublePendulum,solStruct_CD_DoublePendulum, ...
            solStruct_seq_DoublePendulum_regularization);
    end

    if enable_plot_convergence
        % Evaluate
        plotConvergence(probStruct_DoublePendulum, ...
            enable_second_sequential,...
            evalSequentialStruct_DoublePendulum_regularization,...
            evalSequentialStruct_DoublePendulum_hessApprox2,...
            solStruct_CD_DoublePendulum)
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
    solStruct_seq_Satellite_regularization = solve_ctrlDesign_seq_bench_fun(probStruct_Satellite,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_Satellite_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_Satellite_regularization                 = evalSequential(solStruct_seq_Satellite_regularization);
    evalSequentialStruct_Satellite_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Satellite_regularization.probNumber      = kk;

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_Satellite_hessApprox2 = solve_ctrlDesign_seq_bench_fun(probStruct_Satellite,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_Satellite_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_Satellite_hessApprox2 = evalSequential(solStruct_seq_Satellite_hessApprox2);
        evalSequentialStruct_Satellite_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_Satellite_hessApprox2.probNumber = kk;
    end

    % solve Satellite example with two states using custom coordinate-descent
    solStruct_CD_Satellite = solve_ctrlDesign_CD_bench_fun(probStruct_Satellite,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_Satellite;

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows,probStruct_Satellite,solStruct_CD_Satellite, ...
            solStruct_seq_Satellite_regularization,solStruct_seq_Satellite_hessApprox2);
    else
        tableRows = appendResultRow(tableRows,probStruct_Satellite,solStruct_CD_Satellite, ...
            solStruct_seq_Satellite_regularization);
    end

    if enable_plot_convergence
        % Evaluate
        plotConvergence(probStruct_Satellite, ...
            enable_second_sequential,...
            evalSequentialStruct_Satellite_regularization,...
            evalSequentialStruct_Satellite_hessApprox2,...
            solStruct_CD_Satellite)
    end
end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Results table (all benchmarks, all runs)

% final results table accumulated over all runs and all problems
resultsTable = cell2table(tableRows,'VariableNames',tableVariableNames);
disp(resultsTable)


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
plotPerformanceProfiles(seqMethodStructs, cdSolverStructs, solverNames, commonOpts.maxSolveTime)
