function [A, x] = conjugate_gradient(n)
    A = zeros(n, n);
    b = ones(n, 1);
    x = b;
    for i=1:n
    if i==1
    A(i, i:i+1) = [i, 1];
    else if i == n
    A(i, i-1:i) = [1 i];
    else
    A(i, (i-1):(i+1)) = [1 i 1];
    end
    end
    kappa = cond(A);
    res_k = b - A *x;
    p_k = res_k;
    
    res_old = res_k' * res_k;
    for i = 1:n
    Ap = A * p_k;
    alpha = res_old / (p_k' * Ap);
    x = x + alpha * p_k;
    res_k = res_k - alpha * Ap;
    res_new = res_k' * res_k;
    if sqrt(res_new) < 1e-10
    break;
    end
    p_k = res_k + (res_new / res_old) * p_k;
    res_old = res_new;
    it(i) = i;
    res_norm(i) = norm(res_k, 2);
    actual_res_norm(i) = norm(b - Ax, 2);
    estimate(i) = 2 * ((sqrt(kappa) - 1) / (sqrt(kappa) + 1))^i;
    end
    plot(it, res_norm, "o", it, actual_res_norm, it ,estimate, "LineWidth", 1.5);
    xlabel("Iteration", "fontsize", 14)
    ylabel("||r_n||, ||b-Ax_n||, Estimate", "fontsize", 14)
    lengend("Residual norm", "Actual residual norm", "Estimate")
    title("Iteration-wise norms and estimates", "fontsize", 14)
end

[A, x] = conjugate_gradient(100); 
