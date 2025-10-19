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
u = eps('double')/2;

res_sylv = zeros(1, n_tests);
% The first dimension (each row) is a different algorithms that we consider:
%   1. sylvester_mprec_reorth (Algorithm 4.1) in TF32
%   2. sylvester_mprec_inv (Algorithm 4.2) in TF32
%   3. sylvester_mprec_reorth (Algorithm 4.1) in 24-bit custom format
%   4. sylvester_mprec_inv (Algorithm 4.2) in 24-bit custom format
res_mprec = zeros(4, n_tests);
iter = zeros(4, n_tests);
cond_number = zeros(1,n_tests);
Xmprec = {};
is_lyap = false(1, n_tests);

for i = 1:n_tests

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
  tol = 1e-12 * max(m,n);
  max_it = 20;

  tic
  [Xmprec{1}, iter(1, i)] = sylvester_mprec_reorth(coeff1, coeff2, rhs, tol, max_it, @convert_to_tf32);
  toc
  tic
  [Xmprec{2}, iter(2, i)] = sylvester_mprec_inv(coeff1, coeff2, rhs, tol, max_it, @convert_to_tf32);
  toc
  tic
  [Xmprec{3}, iter(3, i)] = sylvester_mprec_reorth(coeff1, coeff2, rhs, tol, max_it, @convert_to_custom_format);
  toc
  tic
  [Xmprec{4}, iter(4, i)] = sylvester_mprec_inv(coeff1, coeff2, rhs, tol, max_it, @convert_to_custom_format);
  toc

  % Condition numbers
  % cond_number(i) = compute_cond_sylv(coeff1, coeff2);
  cond_number(i) = compute_cond_sylv1(coeff1', coeff2');

  % Print results to screen
  res_sylv(i) = norm(rhs - coeff1*Xsylv - Xsylv*coeff2, 2) /...
      (norm(rhs, 2) + norm(Xsylv, 2)*(norm(coeff1, 2)+norm(coeff2, 2)));

  for j = 1 : 4
    res_mprec(j, i) = norm(rhs - coeff1*Xmprec{j} - Xmprec{j}*coeff2, 2) /...
        (norm(rhs, 2) + norm(Xmprec{j}, 2)*(norm(coeff1, 2)+norm(coeff2, 2)));
  end

  fprintf('sylvester() has residual                                %.2e\n', res_sylv(i));
  fprintf('sylvester_mprec_reorth() in TF32 has residual           %.2e\n', res_mprec(1, i));
  fprintf('sylvester_mprec_inv() in TF32 has residual              %.2e\n', res_mprec(2, i));
  fprintf('sylvester_mprec_reorth() in custom format has residual  %.2e\n', res_mprec(1, i));
  fprintf('sylvester_mprec_inv() in custom format has residual     %.2e\n', res_mprec(2, i));
  fprintf('The condition number of the equation is                 %.2e\n', cond_number(i));
end

%% Plot results.
close

% Sort by condition number.
[~, indices] = sort(cond_number, 'descend');

sylvester_strings = 'bx';
mprec_strings = {'vm', '^g', '>k', '<b'};

for k = [1, 2]
  subplot(2,2,(k-1)*2+1)
  semilogy(1:n_tests, res_sylv(indices), sylvester_strings);
  hold on
  for j = 1:2
    semilogy(1:n_tests, res_mprec((k-1)*2+j, indices), mprec_strings{(k-1)*2+j});
  end
  semilogy(1:n_tests, eps() * cond_number(indices), '-');
  hold off
  axis([0, n_tests+1, 1e-24, 1e-10]);
  legend('sylvester',...
         'sylvester\_mprec\_reorth',...
         'sylvester\_mprec\_inv',...
         'condu');
  legend('Location','northeast');
  title('Residual')
  hold off

  subplot(2,2,2*k)

  hold on
  for j = 1:2
    plot(indices, iter((k-1)*2+j, indices), mprec_strings{(k-1)*2+j});
  end
  hold off
  axis([0, n_tests+1, 0, 20]);
  title ('Number of iterations')
end

%% Save results to files.
filename_suffix = {'tf32', '24bit'};
for k = [1, 2]
  outfilename_sylv = sprintf('%s/%s_%s.dat', datfolder, 'test_mixedprecision_sylv', filename_suffix{k});
  outfile_sylv = fopen(outfilename_sylv, 'w');
  id_sylv = 1;

  outfilename_lyap = sprintf('%s/%s_%s.dat', datfolder, 'test_mixedprecision_lyap', filename_suffix{k});
  outfile_lyap = fopen(outfilename_lyap, 'w');
  id_lyap = 1;

  header = ['id    condu res_sylv  r_or i_or  r_in i_in\n'];
  fprintf(outfile_sylv, header);
  fprintf(outfile_lyap, header);
  for i = indices
    if is_lyap(i)
      outfile = outfile_lyap;
      id = id_lyap;
      id_lyap = id_lyap + 1;
    else
      outfile = outfile_sylv;
      id = id_sylv;
      id_sylv = id_sylv + 1;
    end
    fprintf(outfile, '%2d    %.3e %.3e  ',...
            id, eps() * cond_number(i), res_sylv(i));
    for j = [1, 2]
      fprintf(outfile, '%.3e %2d  ',...
              res_mprec((k-1)*2+j, i), iter((k-1)*2+j, i));
    end
    fprintf(outfile, '\n');

  end
  fclose(outfile_sylv);
  fclose(outfile_lyap);
end

function X = convert_to_tf32(Y)
  fpopts.format = 't';
  fpopts.explim = true;
  fpopts.round = 1;
  X = cpfloat(Y, fpopts);
end

function X = convert_to_custom_format(Y)
  fpopts.format = 'custom';
  fpopts.params = [16, -126, 127];
  fpopts.explim = true;
  fpopts.round = 1;
  X = cpfloat(Y, fpopts);
end
