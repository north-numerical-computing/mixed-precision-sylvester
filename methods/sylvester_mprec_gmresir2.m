function [X, num_iter] = sylvester_mprec_gmresir2(A, B, C, ug, iter_max, gtol, reduce_precision)
  %SYLVESTER_MPREC_GMRESIR2    GMRES-based IR for the Sylvester equation.
  %   X = SYLVESTER_MPREC_GMRESIR2(A,B,C,UG,ITER_MAX,GTOL,REDUCE_PRECISION)
  %     solves AX + XB = C using GMRES-based iterative refinement with at most
  %     ITER_MAX refinement steps and GMRES convergence tolerance GTOL.
 %
  %     The global precision is UG, and the function REDUCE_PRECISION casts to the
  %     low precision of interest. The preconditioners are constructed from the
  %     Schur factors of the coefficient matrices A and B.
  %
  %     The residual is computed in double precision.

  if isequal(ug,'uh')
    ug = 'double';
  elseif isequal(ug,'ul')
    ug = 'single';
  else
    assert((isequal(ug,'double')) || (isequal(ug,'single')))
  end

  [dim1, dim2] = size(C);
  n = dim1 * dim2;

  % Default parameters if not specified.
  if nargin < 5
    iter_max = 20;
  end

  if nargin < 6
    gtol = 1e-10 * max(dim1,dim2);
  end

  if (nargin < 7)
    reduce_precision = @(x)single(x);
  end

  % Compute Schur factorizations in low precision.
  [UA, TA] = schur(A);
  UA = reduce_precision(UA);
  TA = reduce_precision(TA);

  [UB, TB] = schur(B);
  UB = reduce_precision(UB);
  TB = reduce_precision(TB);

  % working precision is double
  u = eps('double')/2;
  A = double(A);
  B = double(B);
  C = double(C);

  % Cast Schur factors to working precision.
  UA = double(UA);
  UB = double(UB);
  TA = double(TA);
  TB = double(TB);

  % Use the zero vector as initial guess.
  x =  zeros(n,1,'double');

  cged = false;
  iter = 0; dx = 0; rd = 0;

  gmresits = []; % Number of iterations in each refinement step.
  gmreserr = []; % Final relative residual in each refinement step.

  while ~cged

    % Compute size of errors, quantities in bounds.
    X = reshape(x, dim1, dim2);
    res = C - A*X - X*B;
    nbe(iter+1) = double(norm(res,'inf') / ((norm(A,'inf')+norm(B,'inf'))*norm(X,'inf')+...
                                            norm(C,'inf')));
    temp = abs(res) ./ ( abs(A)*abs(X)+ abs(X)*abs(B) + abs(C)) ;
    temp(isnan(temp)) = 0; % Set 0/0 to 0.
    cbe(iter+1) = max(temp(:));

    iter = iter + 1;
    if iter > iter_max, break, end

    % Check convergence.
    if max([nbe(iter) cbe(iter)]) <= u, break, end

    % Compute residual.
    rd = C - A*X - X*B;

    % Scale residual.
    norm_rd = norm(rd,inf);
    rd1 = rd/norm_rd;

    % Call one step of GMRES to solve for correction term.
    X0 = zeros(dim1, dim2);
    restart = n;
    if n > 10000,
      restart = 1000; % restart is on if problem size too big
    end
    [d, err, its, ~] = schur_precond_sylv_gmres( A, B, rd1, X0, UA, UB, TA, TB, ug, restart, 1, gtol);
    d = d(:);

    % Record number of iterations GMRES took.
    gmresits = [gmresits,its];

    % Record final relative (preconditioned) residual norm in GMRES.
    gmreserr = [gmreserr,err(end)];

    % Record relative (preconditioned) residual norm in each iteration of
    % GMRES (so we can look at convergence trajectories if need be).
    gmreserrvec{iter} = err;

    xold = x;

    % Update solution.
    x = x + norm_rd*double(d);

    dx = norm(x-xold,'inf')/norm(x,'inf');

    % Check if dx contains infs, nans, or is 0.
    if dx == Inf || isnan(double(dx))
      break;
    end
  end
  X = reshape(x, dim1, dim2);
  num_iter = length(gmresits);
end
