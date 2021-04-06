function y = fit_curve(A, X)
    fprintf('A: %.2f', A);
    fprintf('\n');
    y = X * A * X';
end