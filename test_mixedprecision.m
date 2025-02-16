compute_bound = false;

% Matrices not featuring E.
% 'Orr-Som' was not inculded as it does not feature a B matrix.
matrices = {
    {'sylvester_equations/bai1', 256, 0.01},... 256*256, Extreme sizes in bai11
    {'sylvester_equations/bai1', 256, 0.1},... 256*256
    {'sylvester_equations/bai1', 256, 1},... 256*256
    {'sylvester_equations/bai2', 256, 1, 1},... 256*256
    {'sylvester_equations/bai2', 256, 1, 10},... 256*256
    {'sylvester_equations/bai2', 256, 1, 100},... 256*256
    {'sylvester_equations/benner04', 50, 50},... 50*50, Extreme sizes in benn04
    {'sylvester_equations/benner04', 1000, 1000},... 1000*1000
    {'sylvester_equations/blw07', 1},... 317*317
    {'sylvester_equations/blw07', 2},... 1182*1182
    {'sylvester_equations/blw07', 3},... 1182*1182
    {'sylvester_equations/blw07', 4},... 900*900
    {'sylvester_equations/blw07', 5},... 216*216
    {'sylvester_equations/blt09'},... 675*675
    {'sylvester_equations/bqq05'},... 500*500
    {'sylvester_equations/ex_rand'},... 500*500
    {'sylvester_equations/filter2D'},... 1668*1668
    {'sylvester_equations/rail_1357'},... 1357*1357
    {'sylvester_equations/sep_conv_diff_eq', 31, 6, 0.1, 0, 0, 0},... 31*31, Params in hure92
    {'sylvester_equations/sep_conv_diff_eq', 31, 2, 0.01, 25, 50, 50},... 31*31
    {'sylvester_equations/sep_conv_diff_eq', 63, 4, 0.001, 50, 100, 50},... 63*63
    {'sylvester_equations/slicot', 'eady'},... 598*598
    {'sylvester_equations/slicot', 'CDplayer'},... 120*120
    {'sylvester_equations/slicot', 'fom'},... 1006*1006
    {'sylvester_equations/slicot', 'random'},... 200*200
    {'sylvester_equations/slicot', 'pde'},... 84*84
    {'sylvester_equations/slicot', 'heat-cont'},... 200*200
    {'sylvester_equations/slicot', 'iss'},... 270*270
    {'sylvester_equations/slicot', 'build'},... 48*48
    {'sylvester_equations/slicot', 'beam'},... 348*348
    {'sylvester_equations/wlm13', 256, 256},... 256*256, Extreme sizes in wlm13
           };

n_tests = length(matrices);

res_sylv = zeros(1, n_tests);
% The first dimension (each row) is a different algorithms that we consider:
%   1. sylvester_mprec_reorth (Algorithm 4.1)
%   2. sylvester_mprec_inv (Algorithm 4.2)
%   3. sylvester_mprec_gmresir2 (Algorithm 5.1 with ug = uh)
%   4. sylvester_mprec_gmresir2 (Algorithm 5.1 with ug = ul)
res_mprec = zeros(4, n_tests);
iter = zeros(4, n_tests);

Xmprec = {};
is_lyap = false(1, n_tests);

for i = 1:n_tests

  % Generate coefficients of matrix equation.
  % rng(7);
  % n = 300;
  % coeff1 = 10*randn(n,n);% + eye(n);
  % coeff2 = 10*randn(n,n);% + eye(n);
  % coeff1 = anymatrix('sylvester_equations/ex_rand');
  % coeff1 = full(anymatrix('sylvester_equations/rail1357'));

  % coeff2 = coeff1';
  % n = size(coeff1, 1);

  % Symmetric case.
  % coeff1 = coeff1 + coeff1';
  % coeff2 = coeff2 + coeff2';

  % Generate right-hand side of equation from the solution.
  % Xsol = 1*randn(n,n);
  % In =  eye(n);
  % rhs = coeff1*Xsol + Xsol*coeff2;
  % NN = norm(Xsol,2);

  fprintf("***");
  fprintf("  %s", matrices{i}{:});
  fprintf("\n");
  [coeff1, coeff2, rhs] = anymatrix(matrices{i}{:});
  if isempty(coeff2)
    coeff2 = coeff1';
    is_lyap(i) = true;
  end
  coeff1 = full(coeff1);
  coeff2 = full(coeff2);
  if isempty(rhs)
    rhs = randn(size(coeff1, 1), size(coeff2, 2));
  else
    rhs = full(rhs);
  end

  %% Run the test.
  tic
  if is_lyap(i)
    Xsylv = lyap(coeff1, -rhs);
  else
    Xsylv = lyap(coeff1, coeff2, -rhs);
  end
  toc
  % tic
  % Xorth = sylvester_32_64_orth(coeff1, coeff2, rhs);
  % toc

  [m, n] = size(rhs);
  tol = 1e-10 * max(m,n);
  max_it = 20;

  tic
  [Xmprec{1}, iter(1, i)] = sylvester_mprec_reorth(coeff1, coeff2, rhs, tol);
  toc
  tic
  [Xmprec{2}, iter(2, i)] = sylvester_mprec_inv(coeff1, coeff2, rhs, tol);
  toc
  % tic
  % [Xmprec{3}, iter(3, i)] = sylvester_mprec_gmresir2(coeff1, coeff2, rhs, 'uh', max_it, tol);
  % toc
  % tic
  % [Xmprec{4}, iter(4, i)] = sylvester_mprec_gmresir2(coeff1, coeff2, rhs, 'ul', max_it, tol);
  % toc

  % Print results to screen
  res_sylv(i) = norm(rhs - coeff1*Xsylv - Xsylv*coeff2, 2) /...
      (norm(rhs, 2) + norm(Xsylv, 2)*(norm(coeff1, 2)+norm(coeff2, 2)));

  for j = 1 : 2
    res_mprec(j, i) = norm(rhs - coeff1*Xmprec{j} - Xmprec{j}*coeff2, 2) /...
        (norm(rhs, 2) + norm(Xmprec{j}, 2)*(norm(coeff1, 2)+norm(coeff2, 2)));
  end

  fprintf('sylvester() has residual                        %.2e\n', res_sylv(i));
  fprintf('sylvester_mprec_reorth() has residual           %.2e\n', res_mprec(1, i));
  fprintf('sylvester_mprec_inv() has residual              %.2e\n', res_mprec(2, i));
  % fprintf('sylvester_mprec_gmresir2_uh() has residual      %.2e\n', res_mprec(3, i));
  % fprintf('sylvester_mprec_gmresir2_ul() has residual      %.2e\n', res_mprec(4, i));

end

%% Plot results.
close

sylvester_strings = 'bx';
mprec_strings = {'vm', '^g', '>k', '<b'};

subplot(2,1,1)
semilogy(1:n_tests, res_sylv(1:n_tests), sylvester_strings);
hold on
for j = 1:4
  semilogy(1:n_tests, res_mprec(j, 1:n_tests), mprec_strings{j});
end
hold off
axis([0, n_tests+1, 1e-24, 1e-10]);
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
axis([0, n_tests+1, 0, 20]);
title ('Number of iterations')

%% Save results to files.
outfilename_sylv = sprintf('%s/%s', datfolder, 'test_mixedprecision_sylv.dat');
outfile_sylv = fopen(outfilename_sylv, 'w');
id_sylv = 1;

outfilename_lyap = sprintf('%s/%s', datfolder, 'test_mixedprecision_lyap.dat');
outfile_lyap = fopen(outfilename_lyap, 'w');
id_lyap = 1;

header = ['id    res_sylv  r_or i_or  r_in i_in ',...
          'r_gmres_uh   i_gmres_uh   r_gmres_ul   i_gmres_ul\n'];
fprintf(outfile_sylv, header);
fprintf(outfile_lyap, header);
for i = 1:n_tests
  if is_lyap(i)
    outfile = outfile_lyap;
    id = id_lyap;
    id_lyap = id_lyap + 1;
  else
    outfile = outfile_sylv;
    id = id_sylv;
    id_sylv = id_sylv + 1;
  end
  fprintf(outfile, '%2d    %.3e   %.3e %2d   %.3e %2d   %.3e %2d   %.3e %2d\n',...
          id, res_sylv(i),...
          res_mprec(1, i), iter(1, i), res_mprec(2, i), iter(2, i),...
          res_mprec(3, i), iter(3, i), res_mprec(4, i), iter(4, i));

end
fclose(outfile_sylv);
fclose(outfile_lyap);

return

%% Plot the eigenvalues.
clf
eigcoeff1 = eig(coeff1);
eigcoeff2 = eig(coeff2);
plot(real(eigcoeff1), imag(eigcoeff1), 'b+');
hold on
plot(real(eigcoeff2), imag(eigcoeff2), 'ro');
