function [X, error, its, flag] = schur_precond_sylv_gmres(A, B, C, X0, UA, UB, TA, TB, ug, restart, max_it, tol)
  %SCHUR_PRECOND_SYLV_GMRES    Left-preconditioned GMRES for Sylvester equation.
  %   SCHUR_PRECOND_SYLV_GMRES(A,B,C,X0,UA,UB,TA,TB,UG,RESTART,MAX_IT,TOL)
  %   solves the Sylvester equation AX + XB = C by using the Generalized Minimal
  %   residual (GMRES) method on the preconditioned linear system
  %   M^{-1}M*vec(X)=M^{-1}*vec(C). The M-by-M matrix A and the N-by-N matrix B
  %   are the square coefficients, the M-by-N matrix C is the right-hand side of
  %   the equation. The M-by-N matrix X0 is an initial guess, UA and UB are the
  %   unitary Schur coefficients of A and B, respectively, and TA and TB are
  %   their triangular Schur factors. The string UG specifies the working
  %   precision, and can be either "uh" (for binary64) or "ul" (for binary32).
  %   RESTART is the number of iteration between restarts, and MAX_IT is the
  %   maximum number of iteration the algorithm will run in the worst case. The
  %   last input argument TOL specifies the tolerance for the stopping
  %   criterion.
  %   [X,ERROR,ITER,FLAG] = SCHUR_PRECOND_SYLV_GMRES(...) returns the M-by-N
  %   solution X. The output parameter ERROR is The norm of the error at the
  %   last iteration, ITER is The number of inner iterations performed, and FLAG
  %   is an integer set to 0 if The algorithm has converged and to 1 otherwise.

  flag = 0;
  its = 0;

  if isequal(ug,'uh')
    ug = 'double';
  elseif isequal(ug,'ul')
    ug = 'single';
  else
    assert((isequal(ug,'double')) || (isequal(ug,'single')))
  end

  [dim1, dim2] = size(C);

  X = double(X0);
  x = X(:);

  % cast Schur factors to double precision
  UA = double(UA);
  UB = double(UB);
  TA = double(TA);
  TB = double(TB);

  % obtain Schur factors of precision ug
  if isequal(ug,'double')
    UA_ug = double(UA);
    UB_ug = double(UB);
    TA_ug = double(TA);
    TB_ug = double(TB);
  else
    UA_ug = single(UA);
    UB_ug = single(UB);
    TA_ug = single(TA);
    TB_ug = single(TB);
  end

  rtmp = C - A*X - X*B; % compute residual in uh (double)
  if isequal(ug,'double')
    rtmp_ug = double(rtmp);
  else
    rtmp_ug = single(rtmp);
  end

  %Apply M^{-1} to the residual at precision ug
  rtmp_ug = UA_ug' * rtmp_ug * UB_ug;
  rtmp_ug = sylvester (TA_ug, TB_ug.', rtmp_ug);
  rtmp_ug = UA_ug * rtmp_ug * UB_ug';
  r = rtmp_ug(:);

  bnrm2 = norm(r);
  if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

  error(1) = norm( r ) / bnrm2;
  if ( error(1) < tol ) return, end

  n = dim1 * dim2;                                  % initialize workspace
  m = restart;
  V(1:n,1:m+1) = zeros(n,m+1, ug); % GMRES run in precision ug
  H(1:m+1,1:m) = zeros(m+1,m, ug);
  cs(1:m) = zeros(m,1, ug);
  sn(1:m) = zeros(m,1, ug);
  e1    = zeros(n,1, ug);
  e1(1) = 1.0;

  for iter = 1:max_it                              % begin iteration
    V(:,1) = r / norm( r );
    s = norm( r )*e1;
    for i = 1:m                     % construct orthonormal basis via GS
      its = its+1;
      vcur = V(:,i);

      vcur = reshape(vcur, dim1, dim2);
      % applies M in precision uh (double)
      vcur = double(vcur);
      vcur = A * vcur + vcur * B;

      %Apply M^{-1} to vcur in precision uh (double)
      vcur = UA' * vcur * UB;
      vcur = sylvester (TA, TB.', vcur);
      vcur = UA * vcur * UB';
      vcur = vcur(:);

      if isequal(ug,'double')
        w = double(vcur);
      else
        w = single(vcur);
      end

      for k = 1:i
        H(k,i)= w'*V(:,k);
        w = w - H(k,i)*V(:,k);
      end
      H(i+1,i) = norm( w );
      V(:,i+1) = w / H(i+1,i);
      for k = 1:i-1                              % apply Givens rotation
        temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
        H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
        H(k,i)   = temp;
      end
      [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); % form i-th rotation matrix
      temp   = cs(i)*s(i);                        % approximate residual norm
      s(i+1) = -sn(i)*s(i);
      s(i)   = temp;
      H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
      H(i+1,i) = 0.0;
      error((iter-1)*m+i+1)  = abs(s(i+1)) / bnrm2;
      if ( error((iter-1)*m+i+1) <= tol )                        % update approximation
        y = H(1:i,1:i) \ s(1:i);                 % and exit
        addvec = V(:,1:i)*y;
        x = x + double(addvec);
        X = reshape(x, dim1, dim2);
        break;
      end
    end

    if ( error(end) <= tol ), break, end
    y = H(1:m,1:m) \ s(1:m);
    addvec = V(:,1:m)*y;
    x = x + double(addvec);                            % update approximation
    X = reshape(x, dim1, dim2);

    % compute residual
    r = C - A*X - X*B;

    % apply preconditioner to the residual
    if isequal(ug,'double')
      r_ug = double(r);
    else
      r_ug = single(r);
    end

    %Apply M^{-1} to the residual at precision ug
    r_ug = UA_ug' * r_ug * UB_ug;
    r_ug = sylvester (TA_ug, TB_ug.', r_ug);
    r_ug = UA_ug * r_ug * UB_ug';
    r = r_ug(:);

    s(i+1) = norm(r);
    error = [error, s(i+1) / bnrm2];
    % check convergence
    if ( error(end) <= tol ), break, end
  end

  if ( error(end) > tol ) flag = 1; end                 % converged

  function [ c, s ] = rotmat( a, b )
  % Compute the Givens rotation matrix parameters for a and b.

    if ( b == 0.0 )
      c = 1.0;
      s = 0.0;
    elseif ( abs(b) > abs(a) )
      temp = a / b;
      s = 1.0 / sqrt( 1.0 + temp^2 );
      c = temp * s;
    else
      temp = b / a;
      c = 1.0 / sqrt( 1.0 + temp^2 );
      s = temp * c;
    end
  end

end
