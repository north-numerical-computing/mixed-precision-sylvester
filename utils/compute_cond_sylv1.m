function res = compute_cond_sylv1(A, B)
  % Compute the 1-norm condition number of the matrix equation AX + XB = C
  % without generating the Kronecker matrix.

  [m, m1] = size(A);
  [n, n1] = size(B);
  assert(m == m1);
  assert(n == n1);
  [UA, TA] = schur(A);
  [UB, TB] = schur(B);

  res = normest1(@Mv, 5) * normest1(@Msolve, 5);

  function res = Mv(flag, x)
    if strcmp(flag, 'dim')
      res = m*n;
    elseif strcmp(flag, 'real')
      res = isreal(A) && isreal(B);
    else
      X = reshape(x, [m, n, length(x(:)) / (m * n)]);
      Y = zeros(size(X));
      for i = 1:size(Y, 3)
        if strcmp(flag, 'notransp')
          Y(:, :, i) = A * X(:, :, i) + X(:, :, i) * B;
        elseif strcmp(flag, 'transp')
          Y(:, :, i) = A' * X(:, :, i) + X(:, :, i) * B';
        end
      end
      res = reshape(Y, size(x));
    end
  end

  function res = Msolve(flag, x)
    if strcmp(flag, 'dim')
      res = m*n;
    elseif strcmp(flag, 'real')
      res = isreal(A) && isreal(B);
    else
      X = reshape(x, [m, n, length(x(:))/(m*n)]);
      Y = zeros(size(X));
      for i = 1 : size(Y, 3)
        if strcmp(flag, 'notransp')
          Y(:, :, i) = UA * lyap(TA, TB, UA' * X(:, :, i) * UB) * UB';
          % Y(:, :, i) = lyap(A, B, X(:, :, i));
        elseif strcmp(flag, 'transp')
          Y(:, :, i) = UA * lyap(TB, TA, UB' * X(:, :, i)' * UA)' * UB';
          % Y(:, :, i) = lyap(A', B', X(:, :, i));
        end
      end
      res = reshape(Y, size(x));
    end
  end

end
