function Q = myorth(A)
  %MYORTH    Orthogonal factor of QR factorization.
  %   MYORTH(A) is the orthogonal factor of the QR decomposition of A
  %   that corresponds to the matrix R with all positive entries along
  %   the diagonal.
  [Q, R] = qr(A);
  Q = Q .* sign(diag(R))';
end
