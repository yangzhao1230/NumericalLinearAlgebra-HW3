function m=q7(k)
  X1 = rand(2^k, 2^k);
  X2 = rand(2^k, 2^k);
  m = matrix_mul(X1, X2);
  disp(m)
end

function m=matrix_mul(X1, X2)
  s = size(X1, 1)
  if s == 1:
    m = X1 * X2;
  else
    sub_size = s/2;
    A = X1(1:sub_size, 1:sub_size);
    B = X1(1:sub_size, sub_size+1:s);
    C = X1(sub_size+1:s, 1:sub_size);
    D = X1(sub_size+1:s, sub_size+1:s);

  E = X2(1:sub_size, 1:sub_size);
  F = X2(1:sub_size, sub_size+1:s);
  G = X2(sub_size+1:s, 1:sub_size);
  H = X2(sub_size+1:s, sub_size+1:s);
  P1 = matrix_mul(A+D, E+H);
  P2 = matrix_mul(C+D, E);
  P3 = matrix_mul(A, F-H);
  P4 = matrix_mul(D, G-E);
  P5 = matrix_mul(A+B, H);
  P6 = matrix_mul(C-A, F+E);
  P7 = matrix_mul(B-D, G+H);
  W = P1+P4-P5+P7;
  X = P3+P5;
  Y = P2+P4;
  Z = P1+P3-P2+P6;
  m = [W, X; Y, Z];
end
