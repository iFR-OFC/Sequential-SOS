% relative performance profile: fraction of problems per solver column
% whose ratio lies below each grid level,
% f_sx(k, s) = #{ p in P : data(p, s) < x_grid(k) } / n_p
% (strict comparison, as in ClarabelBenchmarks/src/performance_profile.jl;
%  Inf entries of failed solves never count)
function f_sx = perf_profile_cdf(data, x_grid)
% logical 3-D array with entry (p, s, k) = (data(p, s) < x_grid(k)),
% vectorized over all solvers via implicit expansion
f_sx = squeeze(sum(data < reshape(x_grid, 1, 1, []), 1) / size(data, 1)).';
end