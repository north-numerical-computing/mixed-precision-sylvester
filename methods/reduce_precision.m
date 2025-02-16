function X = reduce_precision(Y)
  %REDUCE_PRECISION    Convert matrix to lower precision.
  %   REDUCE_PRECISION(Y) has the same entries as X, converted to lower
  %   binary32 using round-to-nearest with ties-to-even.
  fpopts.format = 's';
  fpopts.explim = true;
  fpopts.round = 1;
  X = cpfloat(Y, fpopts);
end
