function [X, iter] = sylvester_mprec_reorth(A, B, C, tol, max_it, reduce_precision)
  % SYLVESTER_MPREC_REORTH Solve sylvester equation using mixed precision.
  %    X = SYLVESTER_MPREC(A,B,C) solves the Sylvester equation A*X + X*B = C
  %    using a mixed-precision altorithm that combines a low-precision Schur
  %    factorization with a high-precision iterative refinement scheme.
  %
  %    SYLVESTER_MPREC_REORTH(A,B,C,TOL) runs the iterative refinement scheme
  %    until the relative magnitude of the update is below TOL. This parameter
  %    defaults to EPS() * MAX(M,N), where EPS() is the machine precision of
  %    binary64 arithmetic and M and N are the order of A and B, respectively.
  %
  %    SYLVESTER_MPREC_REORTH(A,B,C,TOL,MAX_IT) runs the iterative refinement
  %    scheme for no more than MAX_IT iterations. This default is 20.
  %
  %    SYLVESTER_MPREC_REORTH(A,B,C,TOL,MAX_IT,REDUCE_PRECISION) uses the
  %    REDUCE_PRECISION function handle to convert from working precision to
  %    low precision. By default, the low precision arithmetic is binary32.
  %
  %    [X,ITER] = SYLVESTER_MPREC_REORTH(...) returns the number of iterations
  %    needed to refine the solution to tolerance TOL.

  [m, n] = size(C);

  if (nargin < 4)
    tol = 1e-10 * max(m,n);
  end

  if (nargin < 5)
    max_it = 20;
  end

  if (nargin < 6)
    reduce_precision = @(x)single(x);
  end

  % Compute Schur factorizations in low precision.
  [UA, TA] = schur(A);
  UA = reduce_precision(UA);
  TA = reduce_precision(TA);

  [UB, TB] = schur(B);
  UB = reduce_precision(UB);
  TB = reduce_precision(TB);

  % Solve Sylvester equation in low precision.
  Y0 = matlab.internal.math.sylvester_tri(TA, TB, UA'*C*UB,...
                                          'I', 'I', 'notransp');

  Y0 = reduce_precision(Y0);

  % Orthonormalize unitary factor of Schur decomposition in high precision.
  UA_orth = myorth(UA);
  UB_orth = myorth(UB);

  % Run iteration in higher precision.
  Y = Y0;
  LA = UA_orth' * A * UA_orth;
  LB = UB_orth' * B * UB_orth;
  F = UA_orth' * C * UB_orth;

  % Iteration.
  % DX is the current update and X is the current solution.
  E = Y;
  normY = norm(Y, 2);
  iter = 0;
  while (norm(E, 2) / normY > tol) && iter < max_it
    rhs = F - LA*Y - Y*LB;
    E = matlab.internal.math.sylvester_tri(TA, TB, rhs,...
                                           'I', 'I', 'notransp');

    Y = Y + E;
    iter = iter + 1;
  end

  % Recover solution.
  X = UA_orth * Y * UB_orth';

end
