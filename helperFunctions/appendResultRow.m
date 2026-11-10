% append one row to the results table for a single run:
% compares coordinate descent (CD) and sequential (SQ) and stores the
% problem sizes of the underlying SDP (size of A and of the Gram matrices).
% The first column of each row is the benchmark number ('Prob. No.', as '#<no>')
%
% If a fifth argument is given (sequential run with the second Hessian
% estimation method, cf. enable_second_sequential in the runner), the row
% is extended by the build time, solve time, iterations, and final cost
% of that run.
%
% If the sixth argument includeSolved is true, a solved column per runner
% is appended (1 if the run was successful/feasible, 0 otherwise) instead
% of the unified return status.
function tableRows = appendResultRow(tableRows,probStruct,solStruct_cd,solStruct_seq,solStruct_seq_hessApprox2,includeSolved)

% optional: sequential run with the second Hessian estimation method
if nargin < 5
    solStruct_seq_hessApprox2 = [];
end

% optional: append a solved flag (1/0) per runner
if nargin < 6
    includeSolved = false;
end

% guard: if the fifth argument is a logical, it is the includeSolved flag
% (no second sequential run was passed)
if islogical(solStruct_seq_hessApprox2)
    includeSolved = solStruct_seq_hessApprox2;
    solStruct_seq_hessApprox2 = [];
end

% optional: coordinate descent is not implemented for every benchmark
% (e.g., CBF-CLF); pass [] as solStruct_cd to skip the CD columns
hasCD = ~isempty(solStruct_cd);

% coordinate descent
if hasCD
    solveTimeCD  = solStruct_cd.stats.totalSolveTime;
    iterationsCD = solStruct_cd.stats.iterations;
    buildTimeCD  = solStruct_cd.stats.totalBuildTime;
    finalCostCD  = solStruct_cd.stats.cost(end);
    timePerIterCD = solStruct_cd.stats.time_per_iteration;
end

% sequential
solveTimeSeq  = solStruct_seq.stats.total_solve_time;
buildTimeSeq  = solStruct_seq.stats.solver_build_time;
iterationsSeq = solStruct_seq.stats.iterations;
finalCostSeq  = solStruct_seq.stats.single_iterations{end}.cost;
timePerIterSeq = solStruct_seq.stats.time_per_iteration;

% size of the underlying SDP matrix A, i.e. number of constraints x number
% of decision variables (integer type, so that the table shows e.g. 221260
% instead of the scientific notation 2.2126e+05)
n_con = int64(solStruct_seq.info.sossol_convex.sdpsol.conic.size_a(1));
n_dec = int64(solStruct_seq.info.sossol_convex.sdpsol.conic.size_a(2));

% size of the Gram matrix, i.e. one entry per SOS constraint (joined into
% one string since the number of constraints varies across problems)
sizeGram = strjoin(string(solStruct_seq.info.sossol_convex.gram.Kc(:)).', ', ');

% append the row for this run (sorted by problem and run order)
% first column: number of the benchmark (row index of the results table),
% written as '#<no>'
probNo = size(tableRows,1) + 1;
row = {sprintf('#%d',probNo), probStruct.Name, n_con, n_dec, sizeGram};
if hasCD
    row = [row, buildTimeCD, buildTimeSeq, solveTimeCD, solveTimeSeq, ...
        iterationsCD, iterationsSeq, timePerIterCD, timePerIterSeq, finalCostCD, finalCostSeq];
else
    row = [row, buildTimeSeq, solveTimeSeq, iterationsSeq, timePerIterSeq, finalCostSeq];
end

% solved flag per runner: 1 if the run was successful/feasible, 0 otherwise
if includeSolved
    solvedSeq = double(strcmp(solStruct_seq.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS'));
    if hasCD
        solvedCD = double(strcmp(solStruct_cd.stats.UNIFIED_RETURN_STATUS,  'SOLVER_RET_SUCCESS'));
        row = [row, solvedCD, solvedSeq];
    else
        row = [row, solvedSeq];
    end
end

% sequential with the second Hessian estimation method (if activated)
if ~isempty(solStruct_seq_hessApprox2)
    solveTimeSeq2  = solStruct_seq_hessApprox2.stats.total_solve_time;
    buildTimeSeq2  = solStruct_seq_hessApprox2.stats.solver_build_time;
    iterationsSeq2 = solStruct_seq_hessApprox2.stats.iterations;
    finalCostSeq2  = solStruct_seq_hessApprox2.stats.single_iterations{end}.cost;
    timePerIterSeq2 = solStruct_seq_hessApprox2.stats.time_per_iteration;

    row = [row, buildTimeSeq2, solveTimeSeq2, iterationsSeq2, timePerIterSeq2, finalCostSeq2];

    if includeSolved
        solvedSeq2 = double(strcmp(solStruct_seq_hessApprox2.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS'));
        row = [row, solvedSeq2];
    end
end

tableRows(end+1, :) = row;
end
