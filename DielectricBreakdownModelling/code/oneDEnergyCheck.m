syms y(x) eps(x);
eps(x) = 1;
dy = diff(y, x);
ode = diff(eps * dy, x) == 0;
cond = [ dy(0) == 2,  dy(1) == 1];
ySol(x) = dsolve(ode, cond);


x0 = 0:0.1:1;
y0 = double(ySol(x0));
hold on;
%plot(x0, y0);

yp = diff(ySol, x);
%yp0 = yp(x0);

areaPerm = int(yp, 0, 1)
xlim([0, 1]);
fplot(ySol(x), [0 3], 'DisplayName','y', 'LineWidth', 4);
fplot(yp(x), [0 3],  'DisplayName', 'dy/dx', 'LineWidth', 4);
title("Solution to Laplace Equation with \epsilon = exp(x)", 'FontSize', 18);
text(0.5, 0.1, strcat('Dirichlet Energy : ', num2str(double(areaPerm))), 'FontSize', 14);
legend