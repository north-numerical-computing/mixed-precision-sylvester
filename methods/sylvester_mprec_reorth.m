function [X, iter] = sylvester_mprec_reorth(A, B, C, tol, max_it)
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
%    scheme for no more than MAX_IT iterations. This parameter defaults to 20.
%
%    [X,ITER] = SYLVESTER_MPREC_REORTH(...) returns the number of iterations
%    needed to refine the solution to tolerance TOL.

  debug = false;

  [m, n] = size(C);

  if (nargin < 4)
    tol = 1e-10 * max(m,n);
  end

  if (nargin < 5)
    max_it = 20;
  end

  % Compute Schur factorizations in low precision.
  [UA, TA] = schur(A);
  UA = reduce_precision(UA);
  TA = reduce_precision(TA);

  [UB, TB] = schur(B);
  UB = reduce_precision(UB);
  TB = reduce_precision(TB);

  if debug
    subplot(2,2,1)
    imagesc(log10(abs(TA))); colorbar

    subplot(2,2,2)
    imagesc(log10(abs(TB))); colorbar

    subplot(2,2,3)
    imagesc(log10(abs(UA'*C*UB))); colorbar
  end

  % Solve Sylvester equation in low precision.
  Y0 = matlab.internal.math.sylvester_tri(TA, TB, UA'*C*UB,...
                                          'I', 'I', 'notransp');
  % max(abs(TA(:)))
  % max(abs(UA(:)))
  % max(abs(TB(:)))
  % max(abs(UB(:)))
  % max(abs(X0(:)))
  % keyboard

  Y0 = reduce_precision(Y0);

  if debug
    % any(~isfinite(UA(:)))
    % any(~isfinite(UB(:)))

    % any(~isfinite(TA(:)))
    % any(~isfinite(TB(:)))

    % any(~isfinite(X0(:)))
    subplot(2,2,4)
    spy(~isfinite(X0))
  end

  % keyboard

  % Orthonormalize unitary factor of Schur decomposition in high precision.
  UA_orth = myorth(UA);
  UB_orth = myorth(UB);

  % norm(UB_orth'*UB_orth - eye(size(UB_orth)),2)
  % norm(UA_orth'*UA_orth - eye(size(UA_orth)),2)


  % Run iteration in higher precision.
  Y = Y0;
  LA = UA_orth' * A * UA_orth;
  LB = UB_orth' * B * UB_orth;
  F = UA_orth' * C * UB_orth;

  % Iteration.
  % DX is the current update and X is the current solution.
  E = Y;
  normY = norm(Y, 'fro');
  iter = 0;
  while (norm(E, 'fro') / normY > tol) && iter < max_it
    rhs = F - LA*Y - Y*LB;
    E = matlab.internal.math.sylvester_tri(TA, TB, rhs,...
                                           'I', 'I', 'notransp');

    % % Function that computes Ax = vec(TA * unvec(x) + unvec(x) * TB).
    % afun = @(X)reshape(...
    %     TA*reshape(X, [m, n]) + reshape(X, [m, n])*TB,...
    %     [m*n,1]);

    % % Function that computes M\x = vec(kron(I, TA) + kron(I, TB) \ unvec(x)).
    % mfun = @(X)reshape(...
    %     matlab.internal.math.sylvester_tri(TA, TB,...
    %                                        reshape(X, [m, n]),...
    %                                        'I', 'I', 'notransp'),...
    %     [m*n,1]);

    % % Compute update using GMRES.
    % [x,flag,relres,itnum,resvec] = gmres(afun, rhs(:), [], 1e-13, 100, mfun);
    % DX = reshape(x, [m, n]);

    if debug
      norm(TA * E + E * TB - rhs, 1) / norm(E, 1)
    end

    Y = Y + E;
    iter = iter + 1;
  end

  if debug  && iter > 30
    keyboard
  end

  % Recover solution.
  X = UA_orth * Y * UB_orth';

end
