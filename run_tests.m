% Set up the environment.
initialization

% Produce Figure 1.
plot_optk

% Produce Figure2 2 and 3.
test_mixedprecision

% Produce Figure 4.
reduce_precision = @(X)single(X);
test_conditioning
