%--------------------------------------------------------------------------
%
%   Short Descirption:  Execute the benchmark for the N-link robot arm or
%   just use the .mat file for reproduction (recommended). Note, running
%   the full benchmark takes several hours depending on the computing
%   hardware. The default value for evaluteOnlyFlag is true. To rerun, set
%   to false.
%
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%--------------------------------------------------------------------------


clc
close all
clear

addpath("robot/")
deg    = 2;
noRuns = 1;
Nmax   = 4;  % n = 2*N 

% set to false to re-run benchmark
evaluteOnlyFlag = false;
% set to true to generate tikz figure using MATLAB2TIKZ
genTikz = false;

if ~evaluteOnlyFlag
startTotal = tic;

% run sequential SOS
[fval_array_Seq, solverTimes_total_Seq,buildTimes_Seq,iteration_array_Seq,solStatus_array_Seq,solverStats_Seq] = ROA_Nlink_sequential(noRuns,deg,Nmax);

% run coordinate descent
[fval_array_CD, solverTimes_total_CD,buildTimes_CD,iteration_array_CD,solverStats_CD] = ROA_Nlink_CD(noRuns,deg,Nmax);

% run sequential SOS bad initial
[fval_array_Seq_badInit, solverTimes_total_Seq_badInit,buildTimes_Seq_badInit,iteration_array_Seq_badInit,solStatus_array_Seq_badInit,solverStats_Seq_badInit] = ROA_Nlink_sequential_badInit(noRuns,deg,Nmax);

totalBenchTime = toc(startTotal);


% post-processing of coordinate-descent, because some runs are infeasible,
% and shall not appear in the computational statistics
% see ROA_Nlink_CD.m; if infeasible cost is set to -1, so we only take the
% once unequal to minus one.

idx_infes_cd = fval_array_CD ~=-1;

end

%% evaluation
% in case only evaluate use the provided matfile; otherwise take the fresh
% results 
if evaluteOnlyFlag
    load Full_WS_withCDGratitude.mat
end


runs = 2:Nmax;

% plot solve times
figure('Name', 'Realtive improvement regarding solve times and iterations')
subplot(211)
hold on
rel_solveTime = (solverTimes_total_CD(idx_infes_cd)'-solverTimes_total_Seq')./solverTimes_total_CD(idx_infes_cd)'*100;
mean_solveTimeRel = mean(rel_solveTime);
yline(mean_solveTimeRel)
hold on
bar(2*runs,rel_solveTime);
% xlabel('Number of states')
legend('Mean value','rel. Improvement','Location','best')
ylabel('Relative improvment in solve time in percent')
xticks((2:Nmax)*2) 

% plot iterations (no subiterations of feasibiltity restoration)
subplot(212)
hold on
rel_iteration = (iteration_array_CD(idx_infes_cd)'-iteration_array_Seq')./iteration_array_CD(idx_infes_cd)'*100;
mean_rel_iteration = mean(rel_iteration);
yline(mean_rel_iteration)
hold on
bar(2*runs,rel_iteration);
xlabel('Number of states')
ylabel('Relative improvment in iterations in percent')
legend('Mean value','rel. Improvement','Location','best')
xticks((2:Nmax)*2) 


if genTikz 
cleanfigure();
matlab2tikz('NlinkRobot.tex','width','\figW','height','\figH')
end


% plot solve times
figure('Name', 'Solve times and iterations')
subplot(211)
semilogy(runs*2,solverTimes_total_Seq,'-o')
hold on
semilogy(runs(idx_infes_cd)*2,solverTimes_total_CD(idx_infes_cd),'-^')
xlabel('Number of states')
ylabel('Solver time [s]')
legend('sequential','coordinate-descent','Location','northwest')
xticks((2:Nmax)*2) 

% plot iterations (no subiterations of feasibiltity restoration)
subplot(212)
 plot(runs*2,iteration_array_Seq,'-o')
hold on
plot(runs(idx_infes_cd)*2,iteration_array_CD(idx_infes_cd),'-^')
xlabel('Number of states')
ylabel('Iterations [-]')
legend('sequential','coordinate-descent','Location','northwest')
xticks((2:Nmax)*2) 





% time per iteration
figure(3)
plot(runs*2,solverTimes_total_Seq./iteration_array_Seq,'-o')
hold on
%plot(runs*2,solverTimes_total_Seq_badInit./iteration_array_Seq_badInit,'-o')
plot(runs(idx_infes_cd)*2,solverTimes_total_CD(idx_infes_cd)./iteration_array_CD(idx_infes_cd),'-^')

xlabel('Number of states')
ylabel('Time/iteration [s]')
legend('sequential','coordinate-descent','Location','northwest')
xticks((2:Nmax)*2) 

% get number of linear constraints and decision variables (sequential)
n_con = nan(length(solverStats_Seq),1);
n_dec = nan(length(solverStats_Seq),1);
for j = 1:length(solverStats_Seq)
   n_con(j) = solverStats_Seq{j}{1}{end}.conic.size_A(1);
   n_dec(j) = solverStats_Seq{j}{1}{end}.conic.size_A(2);
end


%% generate tables
nStates = 2*runs';
iteration  = iteration_array_Seq';
times = solverTimes_total_Seq';

Table_seq = table(nStates,n_con,n_dec,iteration,times);


iteration  = iteration_array_CD;
times = solverTimes_total_CD;

Table_cd = table(nStates,n_con,n_dec,iteration,times);


iteration  = iteration_array_Seq_badInit';
times = solverTimes_total_Seq_badInit';

Table_seq_bad = table(nStates,n_con,n_dec,iteration,times);


disp(Table_seq)

disp(Table_cd)

disp(Table_seq_bad)