function [] = plotConvergence(problemStruct, ...
                              enable_second_sequential,...
                              evalSequentialStruct,...
                              evalSequentialStruct_hessApprox2,...
                              solStruct_CD)


% Plot convergence diagnostics
figure('Name',[problemStruct.Name ' Benchmark run ' num2str(evalSequentialStruct.probNumber) ])
tiledlayout(1,2)
nexttile
semilogy(abs(evalSequentialStruct.cost(2:end)-evalSequentialStruct.cost(1:end-1)),'LineWidth',1.2), grid on, title('Cost'), hold on,

if enable_second_sequential
    semilogy(abs(evalSequentialStruct_hessApprox2.cost(2:end)-evalSequentialStruct_hessApprox2.cost(1:end-1)),'LineWidth',1.2), hold on,
end

semilogy(abs(solStruct_CD.stats.cost(2:end)-solStruct_CD.stats.cost(1:end-1)),'LineWidth',1.2), hold on,
xlabel('Iterations')
ylabel('$ |f(\xi_{k+1}) - f(\xi_{k}) |$')
if enable_second_sequential
    legend(['Sequential' ' (' evalSequentialStruct.hess_est_method ')'],['Sequential' ' (' evalSequentialStruct_hessApprox2.hess_est_method ')'],'Coordinate-descent')
else
    legend(['Sequential' ' (' evalSequentialStruct.hess_est_method ')'],'Coordinate-descent')
end

% constraint violation of sequential SOS
% ------------------------------------------------------------------------
% Note: If we are in the SOS cone, constraint violation is negative if
% signed distance is used. In that case, constraint violation is set to
% zero. For plotting on logarithmic scale zero values are set to eps!
% ------------------------------------------------------------------------
nexttile

if evalSequentialStruct.iterations >= evalSequentialStruct_hessApprox2.iterations
    maxIter = evalSequentialStruct.iterations;
else
    maxIter = evalSequentialStruct_hessApprox2.iterations;
end
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
legend({['Sequential' ' (' evalSequentialStruct.hess_est_method ')'],['Sequential'  ' (' evalSequentialStruct_hessApprox2.hess_est_method ')'],'$\epsilon_\mathrm{con}$','eps'})
else
    legend({['Sequential' '(' evalSequentialStruct.hess_est_method ')'],'$\epsilon_\mathrm{con}$','eps'})
end



ax = gca;
ax.YScale = 'log';

end