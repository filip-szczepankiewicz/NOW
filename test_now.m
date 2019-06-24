clear

p = optimizationProblem();

p.N = 100;
p.useMaxNorm = 0;
p.doMaxwellComp = 1;

p = optimizationProblem(p);

mt = tic();
[r, p] = NOW_MULTISCALE(p);
toc(mt)


now_plot_all(r)