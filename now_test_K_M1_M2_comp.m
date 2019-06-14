clear

p = optimizationProblem();

p.N = 50;
p.redoIfFailed = 0;
p.useMaxNorm = 1;

p.durationSecondPartRequested = 24;

p.FlowIndex = 1;
p.AccIndex  = 1;

switch 2
    case 1 %STE
        p.MaxwellIndex = 10^10;
        p.KmatrixIndex = 1000;
        
        p.targetTensor = [1 0 0; 0 1 0; 0 0 1];
    case 2 %LTE
        p.MaxwellIndex = 100;
        p.KmatrixIndex = 10^10;
        
        p.targetTensor = [1 0 0; 0 0 0; 0 0 0];
end

p = optimizationProblem(p);

now_print_requested_and_real_times(p)

[r, p] = NOW_MULTISCALE(p, [20 50], [3 1]);

zind = (diag(p.targetTensor) == 0)';
r.g(:,zind) = 0;
r.gwf(:,zind) = 0;

%%
figure (1)
now_plot_all(r)

figure(2)
f = msf_const_gamma * cumsum(r.g / 1000 .* ( [1 1 1]' * (linspace(0,p.N, p.N+1) * r.dt) )', 1) * r.dt;
a = msf_const_gamma * cumsum(r.g / 1000 .* ( [1 1 1]' * (linspace(0,p.N, p.N+1) * r.dt).^2 )', 1) * r.dt;


subplot(2,1,1)
plot(f)

subplot(2,1,2)
plot(a)