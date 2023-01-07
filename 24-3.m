A = randn(10) - eye(10)
  
t = linspace(0,20)
for i=1:length(t)
    e = expm(t(i)*A)
    y(i) = norm(e)
end
plot(log(t),y)
  
A = randn(10) - eye(10)
t = linspace(0,20)
spectralr = max(abs(eig(A)))
for i=1:length(t)
    e = expm(t(i)*A)
    y(i) = norm(e)
    y1(i) = exp(spectralr*t(i))
end
plot(log(t),y)
hold on;
plot(log(t),y1)
xlabel('log(t)')
ylabel('e(tA),e(a(A)t)')
legend('expm(tA)','exp(a(A)*t)')
