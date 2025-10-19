warning off
compute_bound = false;

cond_magnitudes = [0:0.5:15];
n_tests = length(cond_magnitudes);
cond_nums = zeros(1, n_tests);

res_sylv = zeros(1, n_tests);
% The first dimension (each row) is a different algorithms that we consider:
%   1. sylvester_mprec_reorth (Algorithm 4.1)
%   2. sylvester_mprec_inv (Algorithm 4.2)
%   3. sylvester_mprec_gmresir2 (Algorithm A.1 with ug = uh)
%   4. sylvester_mprec_gmresir2 (Algorithm A.1 with ug = ul)
res_mprec = zeros(4, n_tests);
iter = zeros(4, n_tests);

Xmprec = {};
is_lyap = false(1, n_tests);

for i = 1:n_tests

  fprintf("***");
  fprintf("  %d", cond_magnitudes(i));
  fprintf("\n");

  m = 10;
  n = 10;

  cond_magnitude = cond_magnitudes(i);

  P = randn(m, m);
  coeff1 = P * diag(logspace(0, cond_magnitude, m)) / P;
  P = rand(n, n);
  coeff2 = P * diag(logspace(0, cond_magnitude, m)) / P;
  rhs = randn(m, n);

  fprintf("generated")

  cond_nums(i) = cond(kron(eye(n), coeff1) + kron(coeff2.', eye(m)));

  %% Run the test.
  tic
  if is_lyap(i)
    Xsylv = lyap(coeff1, -rhs);
  else
    Xsylv = lyap(coeff1, coeff2, -rhs);
  end
  toc

  tol = 1e-12 * max(m,n);
  max_it = 20;

  tic
  [Xmprec{1}, iter(1, i)] = sylvester_mprec_reorth(coeff1, coeff2, rhs, tol, max_it, reduce_precision);
  toc
  tic
  [Xmprec{2}, iter(2, i)] = sylvester_mprec_inv(coeff1, coeff2, rhs, tol, max_it, reduce_precision);
  toc
  tic
  [Xmprec{3}, iter(3, i)] = sylvester_mprec_gmresir2(coeff1, coeff2, rhs, 'uh', max_it, tol, reduce_precision);
  toc
  tic
  [Xmprec{4}, iter(4, i)] = sylvester_mprec_gmresir2(coeff1, coeff2, rhs, 'ul', max_it, tol, reduce_precision);
  toc

  % Print results to screen
  res_sylv(i) = norm(rhs - coeff1*Xsylv - Xsylv*coeff2, 2) /...
      (norm(rhs, 2) + norm(Xsylv, 2)*(norm(coeff1, 2)+norm(coeff2, 2)));

  for j = 1 : 4
    res_mprec(j, i) = norm(rhs - coeff1*Xmprec{j} - Xmprec{j}*coeff2, 2) /...
        (norm(rhs, 2) + norm(Xmprec{j}, 2)*(norm(coeff1, 2)+norm(coeff2, 2)));
  end

  fprintf('sylvester() has residual                        %.2e\n', res_sylv(i));
  fprintf('sylvester_mprec_reorth() has residual           %.2e\n', res_mprec(1, i));
  fprintf('sylvester_mprec_inv() has residual              %.2e\n', res_mprec(2, i));
  fprintf('sylvester_mprec_gmresir2_uh() has residual      %.2e\n', res_mprec(3, i));
  fprintf('sylvester_mprec_gmresir2_ul() has residual      %.2e\n', res_mprec(4, i));

end

%% Plot results.
close

sylvester_strings = 'bx';
mprec_strings = {'vm', '^g', '>k', '<b'};
condu_string = 'b';

subplot(2,1,1)
loglog(10 .^ cond_magnitudes(1:n_tests), res_sylv(1:n_tests), sylvester_strings);
hold on
for j = 1:4
  loglog(10 .^ cond_magnitudes(1:n_tests), res_mprec(j, 1:n_tests), mprec_strings{j});
end
semilogy(10 .^ cond_magnitudes(1:n_tests), cond_nums(1:n_tests) * eps() / 2, condu_string);
hold off
axis([0, 10 .^ cond_magnitudes(n_tests), 1e-20, 1e-0]);
legend('sylvester',...
       'sylvester\_mprec\_reorth',...
       'sylvester\_mprec\_inv',...
       'sylvester\_mprec\_gmresir2\_uh',...
       'sylvester\_mprec\_gmresir2\_ul');
legend('Location','northeastoutside');
title('Residual')
hold off

subplot(2,1,2)

hold on
for j = 1:4
  plot(1:n_tests, iter(j, 1:n_tests), mprec_strings{j});
end
hold off
axis([0, n_tests+1, 0, max_it]);
title ('Number of iterations')

%% Save results to files.
outfilename = sprintf('%s/%s', datfolder, 'test_conditioning.dat');
outfile = fopen(outfilename, 'w');

header = ['t condu    res_sylv  r_or i_or  r_in i_in ',...
          'r_gmres_uh   i_gmres_uh   r_gmres_ul   i_gmres_ul\n'];
fprintf(outfile, header);
for i = 1:n_tests
  fprintf(outfile, '%.3e   %.3e   %.3e   %.3e %2d   %.3e %2d   %.3e %2d   %.3e %2d\n',...
          cond_magnitudes(i), cond_nums(i) * eps() / 2, res_sylv(i),...
          res_mprec(1, i), iter(1, i), res_mprec(2, i), iter(2, i),...
          res_mprec(3, i), iter(3, i), res_mprec(4, i), iter(4, i));

end
fclose(outfile);

return

%% Plot the eigenvalues.
clf
eigcoeff1 = eig(coeff1);
eigcoeff2 = eig(coeff2);
plot(real(eigcoeff1), imag(eigcoeff1), 'b+');
hold on
plot(real(eigcoeff2), imag(eigcoeff2), 'ro');
