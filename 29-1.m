%% (a)

function T=tridiag(A)
    A = hilb(4)
    T = hess(A)
end

%% (b)

function Tnew=qralg(T)
m = 4
Ttemp=T
while abs(Ttemp(m, m-1)) > 1e-12
[q, r] = qr(Ttemp)
Ttemp = r * q
end
Tnew = Ttemp
end

%% (d)

function mu=WilkinsonShift(a,b,c)
delta = (a - c) / 2
mu = c - sign(delta) * b^2 / (abs(delta) + sqrt(delta^2 + b^2))
end
function Tnew = qralg(T)
m=4
Ttemp = T
I = eye(m)
while abs(Ttemp(m, m-1)) > 1e-12
mu = WilkinsonShift(Ttemp(m-1, m-1), Ttemp(m,m), Ttemp(m-1, m))
[q, r] = qr(Ttemp - mu * I)
Ttemp = r * q + mu * I
4
end
Tnew = Ttemp
end
