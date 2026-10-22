% append one row to the results table for a single run:
% compares coordinate descent (CD) and sequential (SQ) and stores the
% problem sizes of the underlying SDP (size of A and of the Gram matrices)
function tableRows = appendResultRow(tableRows,probStruct,solStruct_cd,solStruct_seq)

% coordinate descent
solveTimeCD  = solStruct_cd.stats.totalSolveTime;
iterationsCD = solStruct_cd.stats.iterations;
buildTimeCD  = solStruct_cd.stats.totalBuildTime;
finalCostCD  = solStruct_cd.stats.cost(end);

% sequential
solveTimeSeq  = solStruct_seq.stats.total_solve_time;
buildTimeSeq  = solStruct_seq.stats.solver_build_time;
iterationsSeq = solStruct_seq.stats.iterations;
finalCostSeq  = solStruct_seq.stats.single_iterations{end}.cost;

% size of the underlying SDP matrix A, i.e. number of constraints x number
% of decision variables
n_con = solStruct_seq.info.sossol_convex.sdpsol.conic.size_a(1);
n_dec = solStruct_seq.info.sossol_convex.sdpsol.conic.size_a(2);

% size of the Gram matrix, i.e. one entry per SOS constraint (joined into
% one string since the number of constraints varies across problems)
sizeGram = strjoin(string(solStruct_seq.info.sossol_convex.gram.Kc(:)).', ', ');

% append the row for this run (sorted by problem and run order)
tableRows(end+1, :) = {probStruct.Name, n_con, n_dec, sizeGram, ...
    buildTimeCD, buildTimeSeq, solveTimeCD, solveTimeSeq, ...
    iterationsCD, iterationsSeq, finalCostCD, finalCostSeq};
end
