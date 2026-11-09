function [] = plotConvergence(problemStruct, ...
                              enable_second_sequential,...
                              evalSequentialStruct,...
                              evalSequentialStruct_hessApprox2,...
                              solStruct_CD)


% Plot convergence diagnostics
figure('Name',[problemStruct.Name ' Benchmark run ' num2str(evalSequentialStruct.probNumber) ])
tiledlayout(1,2)
nexttile
semilogy(0:evalSequentialStruct.iterations,evalSequentialStruct.scaled_d_ldx,'-*','LineWidth',1.2), grid on, title('Cost'), hold on,

if enable_second_sequential
    semilogy(0:evalSequentialStruct_hessApprox2.iterations,evalSequentialStruct_hessApprox2.scaled_d_ldx,'-o','LineWidth',1.2), hold on,
end

% coordinate descent is optional (not implemented for every benchmark,
% e.g. CBF-CLF); pass [] as solStruct_CD to skip its curve
if ~isempty(solStruct_CD)
    semilogy(0:solStruct_CD.stats.iterations,solStruct_CD.stats.nabla_L_scaled,'-+','LineWidth',1.2), hold on,
end
xlabel('Iterations')
ylabel(['$ \frac{\left\| \nabla_{L}(\xi^k,\lambda^k) \right\|_\infty \, \left\| \omega^k \right\|_\infty }' ...
        '{\max\left(1,\left|f(\xi^k)\right|\right) + \left\|\langle \lambda^k, g(\xi^k)\rangle\right\|_\infty} $'], ...
       'Interpreter','latex')

legendLabels = {['Sequential' ' (' evalSequentialStruct.hess_est_method ')']};
if enable_second_sequential
    legendLabels{end+1} = ['Sequential' ' (' evalSequentialStruct_hessApprox2.hess_est_method ')'];
end
if ~isempty(solStruct_CD)
    legendLabels{end+1} = 'Coordinate-descent';
end
legend(legendLabels{:})
if enable_second_sequential
    maxIter = max([evalSequentialStruct.iterations, evalSequentialStruct_hessApprox2.iterations, solStruct_CD.stats.iterations]);
else
    maxIter = max(evalSequentialStruct.iterations,solStruct_CD.stats.iterations);
end
plot([0 maxIter],[evalSequentialStruct.options.tolerance_opt evalSequentialStruct.options.tolerance_opt],'k--','LineWidth',1.2)

% constraint violation of sequential SOS
% ------------------------------------------------------------------------
% Note: If we are in the SOS cone, constraint violation is negative if
% signed distance is used. In that case, constraint violation is set to
% zero. For plotting on logarithmic scale zero values are set to eps!
% ------------------------------------------------------------------------
nexttile

% if enable_second_sequential
%     maxIter = max(evalSequentialStruct.iterations, evalSequentialStruct_hessApprox2.iterations);
% else
%     maxIter = evalSequentialStruct.iterations;
% end
% plot tolerance for constraint satisfaction
plot([0 maxIter],[evalSequentialStruct.options.tolerance_con evalSequentialStruct.options.tolerance_con],'k--','LineWidth',1.2)
hold on
% plot eps line, i.e., our "zero"
plot([0 maxIter],[eps eps],'k-.','LineWidth',1.2)


semilogy(evalSequentialStruct.maxConVio,'LineWidth',1.2), grid on, title('Maximum constraint violation')


if enable_second_sequential
    semilogy(evalSequentialStruct_hessApprox2.maxConVio,'LineWidth',1.2), hold on
end


xlabel('Iterations')
ylabel('$\min(\mathrm{eps},\theta(\xi_{k}))$')

if enable_second_sequential
legend({'$\epsilon_\mathrm{con}$','eps',['Sequential' ' (' evalSequentialStruct.hess_est_method ')'],['Sequential'  ' (' evalSequentialStruct_hessApprox2.hess_est_method ')']})
else
    legend({'$\epsilon_\mathrm{con}$','eps',['Sequential' '(' evalSequentialStruct.hess_est_method ')']})
end



ax = gca;
ax.YScale = 'log';

end