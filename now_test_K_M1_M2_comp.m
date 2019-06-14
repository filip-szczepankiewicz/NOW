clear

p = optimizationProblem();

p.N = 50;
p.redoIfFailed = 0;
p.useMaxNorm = 1;

p.sMax = 80;

p.durationFirstPartRequested = 31-2;
p.durationSecondPartRequested = 31-6-2;
p.durationZeroGradientRequested = 8;

p.FlowIndex = 1;
p.AccIndex  = 1/1000;


% % Very small values
% p.FlowIndex = .01;
% p.AccIndex  = 1/10000;

% stable max values
% p.AccIndex  = 10^2;
% p.FlowIndex = 10^6;

switch 0
    case 1 %STE Mmatrix
        p.MaxwellIndex = 100;
        p.KmatrixIndex = 10^10;
        
        p.targetTensor = [1 0 0; 0 1 0; 0 0 1];
        
    case 1 %STE Kmatrix
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

%%

[r, p] = NOW_MULTISCALE(p, [20 50], [7 1]);

zind = (diag(p.targetTensor) == 0)';
r.g(:,zind) = 0;
r.gwf(:,zind) = 0;

%%
figure (1)
now_plot_all(r)

save_current_fig_to_file('NOW_STE', '', [6 6], 300)

figure(2)
clf
f = msf_const_gamma * cumsum(r.g / 1000 .* ( [1 1 1]' * (linspace(0,p.N, p.N+1) * r.dt) )', 1) * r.dt;
a = msf_const_gamma * cumsum(r.g / 1000 .* ( [1 1 1]' * (linspace(0,p.N, p.N+1) * r.dt).^2 )', 1) * r.dt;


subplot(2,1,1); hold on
plot(f*0, 'k--')
plot(f)
axis tight
ylabel('Flow')

subplot(2,1,2); hold on
plot(a*0, 'k--')
plot(a)
axis tight
ylabel('Acceleration')


P = gwf_to_pars(r.gwf, r.rf, r.dt);
k0 = P.k0
k1 = P.k1
k2 = P.k2