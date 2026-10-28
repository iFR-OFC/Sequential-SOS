%% ***********************************************************************
%
%  Region-of-attraction Estimation Benchmark Test Runner
%
%  Benchmark 1: Van-der-Pol Oscillator
%  Benchmark 2: Multi-Input-Multi-Output
%  Benchmark 3: Double Pendulum
%  Benchmark 4: Aircraft
%  Benchmark 5: Satellite
%  Benchmark 6: F-18
%  Benchmark 7: Network
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

maxNumCompThreads(1);

% add the problem specific folders and the shared helper functions
runnerDir = fileparts(mfilename('fullpath'));
addPaths  = { ...
    fullfile(runnerDir,'VDP'), ...
    fullfile(runnerDir,'MultiInput'), ...
    fullfile(runnerDir,'DoublePendulum'), ...
    fullfile(runnerDir,'Aircraft'), ...
    fullfile(runnerDir,'LotkaVoltera'), ...
    fullfile(runnerDir,'Satellite_Att'), ...
    fullfile(runnerDir,'F18'), ...
    fullfile(runnerDir,'Network'), ...
    fullfile(runnerDir,'..','helperFunctions')};
for ii = 1:numel(addPaths)
    addpath(addPaths{ii})
end

% enable plotting of convergence for each run
enable_plot_convergence  = true;
enable_second_sequential = true;

% Hessian estimation methods of the sequential solver to compare
hess_est_method_1 = 'regularization'; % our default (see also sequentialCommon.m)
hess_est_method_2 = 'levenberg';           % 'BFGS', 'levenberg', 'mirroring'

% common options shared by the methods
commonOpts.maxiter       = 350;
commonOpts.tolerance_opt = 1e-4;
commonOpts.maxSolveTime  = 3600;
commonOpts.solver        = 'mosek';
commonOpts.miniter       = 5; % minimum iterations before an aborted CD run still counts as successful

% MOSEK settings


commonOpts.sdp_solver_opts.mosek_param = struct( ...
    'MSK_IPAR_NUM_THREADS',            1);

commonOpts.sdp_solver_opts.mosek_echo = 0;

% variable names of the results table (one row per run)
tableVariableNames = {'Problem Name','no. con. (SDP)','no. dec. (SDP)','Size Gram matrices', ...
    'Build Time (CD)','Build Time (SQ)','Solve Time (CD)','Solve Time (SQ)', ...
    'Iterations (CD)','Iterations (SQ)','Final Cost (CD)','Final Cost (SQ)', ...
    'Solved (CD)','Solved (SQ)'};

% extend the results table by the sequential run with the second Hessian
% estimation method (cf. enable_second_sequential)
if enable_second_sequential
    tableVariableNames = [tableVariableNames, ...
        ['Build Time (SQ ' hess_est_method_2 ')'], ...
        ['Solve Time (SQ ' hess_est_method_2 ')'], ...
        ['Iterations (SQ ' hess_est_method_2 ')'], ...
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
degV = {2,2:4,1:4,2:6};
degS = {2,2,2:4,2:4};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_VDP = generateVDP_probStruct(degV{kk},degS{kk});
    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_VDP.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_VDP(probStruct_VDP);

    % solve VDP example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_VDP_regularization = solve_ROA_seq_bench_fun(probStruct_VDP,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_VDP_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_VDP_regularization                 = evalSequential(solStruct_seq_VDP_regularization);
    evalSequentialStruct_VDP_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_VDP_regularization.probNumber      = kk;
    

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_VDP_hessApprox2 = solve_ROA_seq_bench_fun(probStruct_VDP,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_VDP_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_VDP_hessApprox2 = evalSequential(solStruct_seq_VDP_hessApprox2);
        evalSequentialStruct_VDP_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_VDP_hessApprox2.probNumber = kk;
    end

    % solve VDP example with two states using custom coordinate-descent
    solStruct_CD_VDP = solve_ROA_CD_bench_fun(probStruct_VDP,commonOpts,x0);

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


%% Multi-Input-Multi-Output

% define the problem instances based on different polynomial degrees
degV = {2, 2:4, 2:4, 2:6};
degS = {2,   2, 2:4, 2:4};

for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_MiMo = generateMiMo_probStruct(degV{kk},degS{kk});
    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_MiMo.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_MiMo(probStruct_MiMo);

    % solve MiMo example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_MiMo_regularization = solve_ROA_seq_bench_fun(probStruct_MiMo,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_MiMo_regularization;

    % extract size and timing statistic from sequential
    evalSequentialStruct_MiMo_regularization = evalSequential(solStruct_seq_MiMo_regularization);
    evalSequentialStruct_MiMo_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_MiMo_regularization.probNumber = kk;

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_MiMo_hessApprox2 = solve_ROA_seq_bench_fun(probStruct_MiMo,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_MiMo_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_MiMo_hessApprox2 = evalSequential(solStruct_seq_MiMo_hessApprox2);
        evalSequentialStruct_MiMo_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_MiMo_hessApprox2.probNumber = kk;
    end

    % solve MiMo example with two states using custom coordinate-descent
    solStruct_CD_MiMo = solve_ROA_CD_bench_fun(probStruct_MiMo,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_MiMo;
    probNames{end+1}       = sprintf('%s (instance %d of %d)', probStruct_MiMo.Name, kk, numel(degV));


    % Generate table

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows,probStruct_MiMo,solStruct_CD_MiMo, ...
            solStruct_seq_MiMo_regularization,solStruct_seq_MiMo_hessApprox2,true);
    else
        tableRows = appendResultRow(tableRows,probStruct_MiMo,solStruct_CD_MiMo, ...
            solStruct_seq_MiMo_regularization,[],true);
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


%% Double Pendulum

% define the problem instances based on different polynomial degrees
degV = {2;2:4;2:4};
degS = {2;2;2:4};

for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_DoublePendulum = generateDoublePendulum_probStruct(degV{kk},degS{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_DoublePendulum.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_DoublePendulum(probStruct_DoublePendulum);

    % solve DoublePendulum example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_DoublePendulum_regularization = solve_ROA_seq_bench_fun(probStruct_DoublePendulum,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_DoublePendulum_regularization;

    % extract size and timing statistic from sequential
    evalSequentialStruct_DoublePendulum_regularization = evalSequential(solStruct_seq_DoublePendulum_regularization);
    evalSequentialStruct_DoublePendulum_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_DoublePendulum_regularization.probNumber = kk;

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_DoublePendulum_hessApprox2 = solve_ROA_seq_bench_fun(probStruct_DoublePendulum,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_DoublePendulum_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_DoublePendulum_hessApprox2 = evalSequential(solStruct_seq_DoublePendulum_hessApprox2);
        evalSequentialStruct_DoublePendulum_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_DoublePendulum_hessApprox2.probNumber = kk;
    end

    % solve DoublePendulum example with two states using custom coordinate-descent
    solStruct_CD_DoublePendulum = solve_ROA_CD_bench_fun(probStruct_DoublePendulum,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_DoublePendulum;
    probNames{end+1}       = sprintf('%s (instance %d of %d)', probStruct_DoublePendulum.Name, kk, numel(degV));

    % Generate table

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows,probStruct_DoublePendulum,solStruct_CD_DoublePendulum, ...
            solStruct_seq_DoublePendulum_regularization,solStruct_seq_DoublePendulum_hessApprox2,true);
    else
        tableRows = appendResultRow(tableRows,probStruct_DoublePendulum,solStruct_CD_DoublePendulum, ...
            solStruct_seq_DoublePendulum_regularization,[],true);
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


%% Aircraft 4D


% define the problem instances based on different polynomial degrees
degV = {2,2:4,1:4};
degS = {2:4,2:4,2:4};

for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Aircraft = generateAircraft_probStruct(degV{kk},degS{kk});
    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_Aircraft.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_Aircraft(probStruct_Aircraft);

    % solve Aircraft example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Aircraft_regularization = solve_ROA_seq_bench_fun(probStruct_Aircraft,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_Aircraft_regularization;

    % extract size and timing statistic from sequential
    evalSequentialStruct_Aircraft_regularization = evalSequential(solStruct_seq_Aircraft_regularization);
    evalSequentialStruct_Aircraft_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Aircraft_regularization.probNumber = kk;

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_Aircraft_hessApprox2 = solve_ROA_seq_bench_fun(probStruct_Aircraft,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_Aircraft_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_Aircraft_hessApprox2 = evalSequential(solStruct_seq_Aircraft_hessApprox2);
        evalSequentialStruct_Aircraft_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_Aircraft_hessApprox2.probNumber = kk;
    end

    % solve Aircraft example with two states using custom coordinate-descent
    solStruct_CD_Aircraft = solve_ROA_CD_bench_fun(probStruct_Aircraft,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_Aircraft;
    probNames{end+1}       = sprintf('%s (instance %d of %d)', probStruct_Aircraft.Name, kk, numel(degV));

    % Generate table

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows,probStruct_Aircraft,solStruct_CD_Aircraft, ...
            solStruct_seq_Aircraft_regularization,solStruct_seq_Aircraft_hessApprox2,true);
    else
        tableRows = appendResultRow(tableRows,probStruct_Aircraft,solStruct_CD_Aircraft, ...
            solStruct_seq_Aircraft_regularization,[],true);
    end

    if enable_plot_convergence
        % Evaluate
        plotConvergence(probStruct_Aircraft, ...
                        enable_second_sequential,...
                        evalSequentialStruct_Aircraft_regularization,...
                        evalSequentialStruct_Aircraft_hessApprox2,...
                        solStruct_CD_Aircraft)
    end

end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))

% %% LotkaVol
% 
% % define the problem instances based on different polynomial degrees
% degV = {2,2:4,1:4};
% degS = {2,2:4,2:4};
% 
% for kk = 1:numel(degV)
% 
%     % generate problem structure (dynamics, polynomial dec. variables)
%     probStruct_LotkaVol = generateLotkaVol_probStruct(degV{kk},degS{kk});
%     fprintf('------------------------------------------------------------------------------------------ \n')
%     fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_LotkaVol.Name, kk, numel(degV))
%     fprintf('------------------------------------------------------------------------------------------\n')
% 
%     % generate an initial guess
%     x0 = generate_initGuess_LotkaVol(probStruct_LotkaVol);
% 
%     % solve LotkaVol example with two states using sequential SOS
%     hess_est_method = hess_est_method_1;
%     solStruct_seq_LotkaVol_regularization = solve_ROA_seq_bench_fun(probStruct_LotkaVol,commonOpts,x0,hess_est_method);
% 
%     seqSolverStructs_regularization{end+1} = solStruct_seq_LotkaVol_regularization;
% 
%     % extract size and timing statistic from sequential
%     evalSequentialStruct_LotkaVol_regularization = evalSequential(solStruct_seq_LotkaVol_regularization);
%     evalSequentialStruct_LotkaVol_regularization.hess_est_method = hess_est_method_1;
%     evalSequentialStruct_LotkaVol_regularization.probNumber = kk;
% 
%     if enable_second_sequential
% 
%         hess_est_method = hess_est_method_2;
%         solStruct_seq_LotkaVol_hessApprox2 = solve_ROA_seq_bench_fun(probStruct_LotkaVol,commonOpts,x0,hess_est_method);
% 
%         seqSolverStructs_hessApprox2{end+1} = solStruct_seq_LotkaVol_hessApprox2;
% 
%         % extract size and timing statistic from sequential
%         evalSequentialStruct_LotkaVol_hessApprox2 = evalSequential(solStruct_seq_LotkaVol_hessApprox2);
%         evalSequentialStruct_LotkaVol_hessApprox2.hess_est_method = hess_est_method_2;
%         evalSequentialStruct_LotkaVol_hessApprox2.probNumber = kk;
%     end
% 
%     % solve LotkaVol example with two states using custom coordinate-descent
%     solStruct_CD_LotkaVol = solve_ROA_CD_bench_fun(probStruct_LotkaVol,commonOpts,x0);
% 
%     cdSolverStructs{end+1} = solStruct_CD_LotkaVol;
% 
%     % Generate table
% 
%     % append the current run to the accumulating results table
%     if enable_second_sequential
%         tableRows = appendResultRow(tableRows,probStruct_LotkaVol,solStruct_CD_LotkaVol, ...
%             solStruct_seq_LotkaVol_regularization,solStruct_seq_LotkaVol_hessApprox2);
%     else
%         tableRows = appendResultRow(tableRows,probStruct_LotkaVol,solStruct_CD_LotkaVol, ...
%             solStruct_seq_LotkaVol_regularization);
%     end
% 
%     if enable_plot_convergence
%         % Evaluate
%         plotConvergence(probStruct_LotkaVol, ...
%             enable_second_sequential,...
%             evalSequentialStruct_LotkaVol_regularization,...
%             evalSequentialStruct_LotkaVol_hessApprox2,...
%             solStruct_CD_LotkaVol)
%     end
% 
% end
% 
% % display the results table accumulated over all runs so far
% disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Satellite

% define the problem instances based on different polynomial degrees
degV = {2,2:4,2:6};
degS = {2,2,2:4};


for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Satellite = generateSatellite_probStruct(degV{kk},degS{kk});
    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_Satellite.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_Satellite(probStruct_Satellite);

    % solve Satellite example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Satellite_regularization = solve_ROA_seq_bench_fun(probStruct_Satellite,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_Satellite_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_Satellite_regularization                 = evalSequential(solStruct_seq_Satellite_regularization);
    evalSequentialStruct_Satellite_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Satellite_regularization.probNumber      = kk;


    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_Satellite_hessApprox2 = solve_ROA_seq_bench_fun(probStruct_Satellite,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_Satellite_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_Satellite_hessApprox2 = evalSequential(solStruct_seq_Satellite_hessApprox2);
        evalSequentialStruct_Satellite_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_Satellite_hessApprox2.probNumber = kk;
    end

    % solve Satellite example with two states using custom coordinate-descent
    solStruct_CD_Satellite = solve_ROA_CD_bench_fun(probStruct_Satellite,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_Satellite;
    probNames{end+1}       = sprintf('%s (instance %d of %d)', probStruct_Satellite.Name, kk, numel(degV));


    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows,probStruct_Satellite,solStruct_CD_Satellite, ...
            solStruct_seq_Satellite_regularization,solStruct_seq_Satellite_hessApprox2,true);
    else
        tableRows = appendResultRow(tableRows,probStruct_Satellite,solStruct_CD_Satellite, ...
            solStruct_seq_Satellite_regularization,[],true);
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



%% F/A-18

% define the problem instances based on different polynomial degrees
degV = {2, 2:4, 2:4};
degS = {2,   2, 2:4};

for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_F18 = generateF18_probStruct(degV{kk},degS{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_F18.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_F18(probStruct_F18);

    % solve F-18 example with seven states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_7D_regularization = solve_ROA_seq_bench_fun(probStruct_F18,commonOpts,x0,hess_est_method);

    % store in a common struct for all sequential instances and problem
    seqSolverStructs_regularization{end+1} = solStruct_seq_7D_regularization;

    % extract size and timing statistic from sequential
    evalSequentialStruct_7D_regularization = evalSequential(solStruct_seq_7D_regularization);
    evalSequentialStruct_7D_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_7D_regularization.probNumber = kk;

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_7D_hessApprox2 = solve_ROA_seq_bench_fun(probStruct_F18,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_7D_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_7D_hessApprox2 = evalSequential(solStruct_seq_7D_hessApprox2);
        evalSequentialStruct_7D_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_7D_hessApprox2.probNumber = kk;
    end

    % solve F-18 example with seven states using custom coordinate-descent
    solStruct_CD_7D = solve_ROA_CD_bench_fun(probStruct_F18,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_7D;
    probNames{end+1}       = sprintf('%s (instance %d of %d)', probStruct_F18.Name, kk, numel(degV));

    % Generate table

    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows,probStruct_F18,solStruct_CD_7D, ...
            solStruct_seq_7D_regularization,solStruct_seq_7D_hessApprox2,true);
    else
        tableRows = appendResultRow(tableRows,probStruct_F18,solStruct_CD_7D, ...
                                    solStruct_seq_7D_regularization,[],true);
    end


    if enable_plot_convergence
        % Evaluate
        plotConvergence(probStruct_F18, ...
                        enable_second_sequential,...
                        evalSequentialStruct_7D_regularization,...
                        evalSequentialStruct_7D_hessApprox2,...
                        solStruct_CD_7D)
    end

end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))

%% 12-State Dynamical System

% generate problem structure (dynamics, polynomial dec. variables)
probStruct_12State= generate12State_probStruct();
fprintf('------------------------------------------------------------------------------------------ \n')
fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_12State.Name, 1, 1)
fprintf('------------------------------------------------------------------------------------------\n')

% generate an initial guess
x0 = generate_initGuess_12State(probStruct_12State);

% solve 12 state dynamical system example with 12 states using sequential SOS
hess_est_method = hess_est_method_1;
solStruct_seq_12D_regularization = solve_ROA_seq_bench_fun(probStruct_12State,commonOpts,x0,hess_est_method);

seqSolverStructs_regularization{end+1} = solStruct_seq_12D_regularization;

% extract size and timing statistic from sequential
evalSequentialStruct_12D_regularization = evalSequential(solStruct_seq_12D_regularization);
evalSequentialStruct_12D_regularization.hess_est_method = hess_est_method_1;
evalSequentialStruct_12D_regularization.probNumber = 1;

if enable_second_sequential

    hess_est_method = hess_est_method_2;
    solStruct_seq_12D_hessApprox2 = solve_ROA_seq_bench_fun(probStruct_12State,commonOpts,x0,hess_est_method);

    seqSolverStructs_hessApprox2{end+1} = solStruct_seq_12D_hessApprox2;

    % extract size and timing statistic from sequential
    evalSequentialStruct_12D_hessApprox2 = evalSequential(solStruct_seq_12D_hessApprox2);
    evalSequentialStruct_12D_hessApprox2.hess_est_method = hess_est_method_2;
    evalSequentialStruct_12D_hessApprox2.probNumber = 1;
end

% solve 12 state dynamical system example with 12 states using custom coordinate-descent
solStruct_CD_12D = solve_ROA_CD_bench_fun(probStruct_12State,commonOpts,x0);

cdSolverStructs{end+1} = solStruct_CD_12D;
probNames{end+1}       = sprintf('%s (instance %d of %d)', probStruct_12State.Name, 1, 1);

% Generate table

% append the current run to the accumulating results table
if enable_second_sequential
    tableRows = appendResultRow(tableRows,probStruct_12State,solStruct_CD_12D, ...
        solStruct_seq_12D_regularization,solStruct_seq_12D_hessApprox2,true);
else
    tableRows = appendResultRow(tableRows,probStruct_12State,solStruct_CD_12D, ...
        solStruct_seq_12D_regularization,[],true);
end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


if enable_plot_convergence
    % Evaluate
    plotConvergence(probStruct_12State, ...
                    enable_second_sequential,...
                    evalSequentialStruct_12D_regularization,...
                    evalSequentialStruct_12D_hessApprox2,...
                    solStruct_CD_12D)
end


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
plotPerformanceProfiles(seqMethodStructs, cdSolverStructs, solverNames, commonOpts.maxSolveTime, probNames)

% remove the added folders from the path again now that the benchmarks
% are done (exactly the folders added at the top of this script)
for ii = 1:numel(addPaths)
    rmpath(addPaths{ii})
end
