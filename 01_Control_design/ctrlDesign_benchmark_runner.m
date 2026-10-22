%% ***********************************************************************
%
%  Region-of-attraction Estimation Benchmark Test Runner
%
%  Benchmark 1: Van-der-Pol Oscillator
%  Benchmark 2: Multi-Input-Multi-Output
%  Benchmark 3: Double Pendulum
%  Benchmark 4: Aircraft
%  Benchmark 5: F-18
%  Benchmark 6: Network
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

% enable plotting of convergence for each run
enable_plot_convergence = true;
enable_second_sequential = true;

% Hessian estimation methods of the sequential solver to compare
hess_est_method_1 = 'regularization'; % our default (see also sequentialCommon.m)
hess_est_method_2 = 'levenberg';      % 'BFGS', 'levenberg', 'mirroring'

% common options shared by the methods
commonOpts.maxiter       = 350;
commonOpts.tolerance_opt = 1e-4;
commonOpts.maxSolveTime  = 3600;

% initialize the cell array that accumulates one row per run of the results table
tableRows = cell(0, 12);

% variable names of the results table (one row per run)
tableVariableNames = {'Problem Name','no. con. (SDP)','no. dec. (SDP)','Size Gram matrices', ...
    'Build Time (CD)','Build Time (SQ)','Solve Time (CD)','Solve Time (SQ)', ...
    'Iterations (CD)','Iterations (SQ)','Final Cost (CD)','Final Cost (SQ)'};


%% VDP

% define the problem instances based on different polynomial degrees
degV = {2,2:4,1:4,2:6};
degK = {1,1,1,1};

% initialize the cell arrays of stored runs (one entry per run)
seqSolverStructs_regularization = {};
seqSolverStructs_hessApprox2 = {};
cdSolverStructs  = {};

for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_VDP = generateVDP_probStruct(degV{kk},degK{kk});
    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_VDP.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_VDP(probStruct_VDP);

    % solve VDP example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_VDP_regularization = solve_ctrlDesign_seq_bench_fun(probStruct_VDP,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_VDP_regularization;

    % extract size and timing statistic from sequential
    evalSequentialStruct_VDP_regularization = evalSequential(solStruct_seq_VDP_regularization);

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_VDP_hessApprox2 = solve_ctrlDesign_seq_bench_fun(probStruct_VDP,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_VDP_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_VDP_hessApprox2 = evalSequential(solStruct_seq_VDP_hessApprox2);
    end

    % solve VDP example with two states using custom coordinate-descent
    solStruct_CD_VDP = solve_ctrlDesign_CD_bench_fun(probStruct_VDP,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_VDP;


    %% Generate table

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows,probStruct_VDP,solStruct_CD_VDP,solStruct_seq_VDP_regularization);

    if enable_plot_convergence
        % Evaluate

        % Plot convergence diagnostics
        figure('Name','VDP Benchmark')
        tiledlayout(1,2)
        nexttile
        hessLegendNames = {hess_est_method_1};
        semilogy(abs(evalSequentialStruct_VDP_regularization.cost(2:end)-evalSequentialStruct_VDP_regularization.cost(1:end-1)),'LineWidth',1.2), grid on, title('Cost'), hold on,
        if enable_second_sequential
            semilogy(abs(evalSequentialStruct_VDP_hessApprox2.cost(2:end)-evalSequentialStruct_VDP_hessApprox2.cost(1:end-1)),'LineWidth',1.2), hold on,
            hessLegendNames{end+1} = hess_est_method_2;
        end
        semilogy(abs(solStruct_CD_VDP.stats.cost(2:end)-solStruct_CD_VDP.stats.cost(1:end-1)),'LineWidth',1.2), hold on,
        xlabel('Iterations')
        ylabel('$ |f(\xi_{k+1}) - f(\xi_{k}) |$')
        hessLegendNames{end+1} = 'Coordinate-descent';
        legend(hessLegendNames)

        % constraint violation of sequential SOS
        % ------------------------------------------------------------------------
        % Note: If we are in the SOS cone, constraint violation is negative if
        % signed distance is used. In that case, constraint violation is set to
        % zero. For plotting on logarithmic scale zero values are set to eps!
        % ------------------------------------------------------------------------
        nexttile, semilogy(evalSequentialStruct_VDP_regularization.maxConVio,'LineWidth',1.2), grid on, title('Maximum constraint violation') , hold on
        plot([0 solStruct_seq_VDP_regularization.stats.iterations],[solStruct_seq_VDP_regularization.stats.options.tolerance_con solStruct_seq_VDP_regularization.stats.options.tolerance_con],'k--','LineWidth',1.2)
        if enable_second_sequential
            semilogy(evalSequentialStruct_VDP_hessApprox2.maxConVio,'LineWidth',1.2), hold on
            plot([0 solStruct_seq_VDP_hessApprox2.stats.iterations],[solStruct_seq_VDP_hessApprox2.stats.options.tolerance_con solStruct_seq_VDP_hessApprox2.stats.options.tolerance_con],'k--','LineWidth',1.2)
        end
        plot([0 solStruct_seq_VDP_regularization.stats.iterations],[eps eps],'k-.','LineWidth',1.2)
        xlabel('Iterations')
        ylabel('$\min(\mathrm{eps},\theta(\xi_{k}))$')
        legend([hessLegendNames(1:end-1), {'$\epsilon_\mathrm{con}$','eps'}])
    end
end