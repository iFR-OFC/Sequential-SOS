% relative performance profile: 

% fraction of problems per solver column
% whose ratio lies below each grid level,
% 

% rho_sx(k, s) = { p in P : data(p, s) < x_grid(k) } / n_p
function rho = perf_profile_cdf(data, x_grid)

% logical array 
rho = squeeze(sum(data < reshape(x_grid, 1, 1, []), 1) / size(data, 1)).';

end