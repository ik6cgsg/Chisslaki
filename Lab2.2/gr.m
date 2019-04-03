fin = fopen('a.out')
n = fscanf(fin, "%i", 1);
x1 = fscanf(fin, "%f", n);
y1_n = fscanf(fin, "%f", n);
y1_d = fscanf(fin, "%f", n);
x2 = fscanf(fin, "%f", n);
y2_n = fscanf(fin, "%f", n);
y2_d = fscanf(fin, "%f", n);
fclose(fin);

f1 = @(x) exp(x);
f2 = @(x) exp(x) ./ (1 - x .^ 2);

plot(x1, y1_n, x1, y1_d, x1, f1(x1), "g");  
grid on;
title("e^x, 3 nodes");
legend("natural spline", "second spline", "original");

figure;

plot(x2, y2_n, x2, y2_d, x2, f2(x2), "g");
title("e^x / (1 - x^2), 3 nodes");
grid on;
legend("natural spline", "second spline", "original");

fin = fopen('b.out')
n = fscanf(fin, "%i", 1);
x1 = fscanf(fin, "%f", n);
y1_n = fscanf(fin, "%f", n);
y1_d = fscanf(fin, "%f", n);
x2 = fscanf(fin, "%f", n);
y2_n = fscanf(fin, "%f", n);
y2_d = fscanf(fin, "%f", n);
fclose(fin);

f1 = @(x) exp(x);
f2 = @(x) exp(x) ./ (1 - x .^ 2);

figure;
plot(x1, y1_n, x1, y1_d, x1, f1(x1), "g");  
grid on
title("e^x, 7 nodes");
legend("natural spline", "second spline", "original");

figure;

plot(x2, y2_n, x2, y2_d, x2, f2(x2), "g");
title("e^x / (1 - x^2), 7 nodes");
grid on;
legend("natural spline", "second spline", "original");

