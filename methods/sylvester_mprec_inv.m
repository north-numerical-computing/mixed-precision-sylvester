function [X, iter] = sylvester_mprec_inv(A, B, C, tol, max_it)
  % SYLVESTER_MPREC_INV    Solve Sylvester equation using mixed precision.
  %    X = SYLVESTER_MPREC_INV(A,B,C) solves the Sylvester equation A*X + X*B= C
  %    using a mixed-precision altorithm that combines a low-precision Schur
  %    factorization with a high-precision iterative refinement scheme.
  %
  %    SYLVESTER_MPREC_INV(A,B,C,TOL) runs the iterative refinement scheme
  %    until the relative magnitude of the update is below TOL. This parameter
  %    defaults to EPS() * MAX(M,N), where EPS() is the machine precision of
  %    binary64 arithmetic and M and N are the order of A and B, respectively.
  %
  %    SYLVESTER_MPREC_INV(A,B,C,TOL,MAX_IT) runs the iterative refinement
  %    scheme for no more than MAX_IT iterations. This parameter defaults to 20.
  %
  %    [X,ITER] = SYLVESTER_MPREC_INV(...) returns the number of iterations
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
  % norm_delta_TA = norm(TA - reduce_precision(TA),2)

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
  % max(abs(Y0(:)))
  % keyboard

  Y0 = reduce_precision(Y0);

  if debug
    % any(~isfinite(UA(:)))
    % any(~isfinite(UB(:)))

    % any(~isfinite(TA(:)))
    % any(~isfinite(TB(:)))

    % any(~isfinite(Y0(:)))
    subplot(2,2,4)
    spy(~isfinite(Y0))
  end

  % Run iteration in higher precision.
 LA = UA' * A / (UA');
 LB = UB \ B * UB;
 F = UA' * C * UB;

 % Iteration.
 % E is the current update and X is the current solution.
 Y = Y0;
 E = Y0;
 normY = norm(Y, 'fro');

 if debug
   norm(LA*Y + Y*LB - UA'*C*UB, 'fro') / norm(UA'*C*UB, 'fro')
 end

 iter = 0;
 while (norm(E, 'fro') / normY > tol) && iter < max_it
   rhs = F - LA*Y - Y*LB;
   E = matlab.internal.math.sylvester_tri(TA, TB, rhs,...
                                          'I', 'I', 'notransp');

   Y = Y + E;

   if debug
     norm(LA*Y + Y*LB - UA'*C*UB, 'fro') / norm(UA'*C*UB, 'fro')
   end

   iter = iter + 1;
 end

   if debug  && iter > 30
     keyboard
   end

   % Recover solution.
   X = UA' \ Y / UB;
   % cond(UA,2)
   % cond(UB,2)
end
