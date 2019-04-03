fin = fopen('f.out')
n_f = fscanf(fin, "%i", 1);
x1_f = fscanf(fin, "%f", n_f);
y1_f = fscanf(fin, "%f", n_f);
x2_f = fscanf(fin, "%f", n_f);
y2_f = fscanf(fin, "%f", n_f);
fclose(fin);

fin = fopen('s.out')
n_s = fscanf(fin, "%i", 1);
x1_s = fscanf(fin, "%f", n_s);
y1_s = fscanf(fin, "%f", n_s);
x2_s = fscanf(fin, "%f", n_s);
y2_s = fscanf(fin, "%f", n_s);
fclose(fin);

fin = fopen('t.out')
n_t = fscanf(fin, "%i", 1);
x1_t = fscanf(fin, "%f", n_t);
y1_t = fscanf(fin, "%f", n_t);
x2_t = fscanf(fin, "%f", n_t);
y2_t = fscanf(fin, "%f", n_t);
fclose(fin);

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
pol = @(x) 10.0018 * x .^ 2 + 7.46844 * x + 0.999996;

% uniform
plot(x1_f, y1_f, 'k', x1_s, y1_s, 'm', x1_f, f1(x1_f), 'g');
title('e^x, uniform');
hold on;
plot(x1, y1_n, x1, y1_d);  


legend('3 nodes', '5 nodes', 'original', 'natural spline', 'second spline');

fin = fopen('b.out')
n = fscanf(fin, "%i", 1);
x1 = fscanf(fin, "%f", n);
y1_n = fscanf(fin, "%f", n);
y1_d = fscanf(fin, "%f", n);
x2 = fscanf(fin, "%f", n);
y2_n = fscanf(fin, "%f", n);
y2_d = fscanf(fin, "%f", n);
fclose(fin);

figure;

plot(x2_f, y2_f, 'k', x2_s, y2_s, 'm', x2_f, f2(x2_f), 'g');
hold on;
plot(x2, y2_n, x2, y2_d);
title('e^x/(1 - x^2), uniform');
grid on;

legend('3 nodes', '5 nodes', 'original',  'natural spline', 'second spline', 'Location', 'NorthWest');

%chebyshev
fin = fopen('ff.out')
n_f = fscanf(fin, "%i", 1);
x1_f = fscanf(fin, "%f", n_f);
y1_f = fscanf(fin, "%f", n_f);
x2_f = fscanf(fin, "%f", n_f);
y2_f = fscanf(fin, "%f", n_f);
fclose(fin);

fin = fopen('ss.out')
n_s = fscanf(fin, "%i", 1);
x1_s = fscanf(fin, "%f", n_s);
y1_s = fscanf(fin, "%f", n_s);
x2_s = fscanf(fin, "%f", n_s);
y2_s = fscanf(fin, "%f", n_s);
fclose(fin);

fin = fopen('tt.out')
n_t = fscanf(fin, "%i", 1);
x1_t = fscanf(fin, "%f", n_t);
y1_t = fscanf(fin, "%f", n_t);
x2_t = fscanf(fin, "%f", n_t);
y2_t = fscanf(fin, "%f", n_t);
fclose(fin);

figure;

plot(x1_f, y1_f, 'k', x1_s, y1_s, 'm', x1_t, y1_t, 'b', x1_f, f1(x1_f), 'g');
title('e^x, chebyshev');
grid on;

legend('3 nodes', '5 nodes', '9 nodes', 'original');

figure;

plot(x2_f, y2_f, 'k', x2_s, y2_s, 'm', x2_t, y2_t, 'b', x2_s, f2(x2_f), 'g');
title('e^x/(1 - x^2), chebyshev');
grid on;

legend('3 nodes', '5 nodes', '9 nodes', 'original');

%user
fin = fopen('fff.out')
n_f = fscanf(fin, "%i", 1);
x1_f = fscanf(fin, "%f", n_f);
y1_f = fscanf(fin, "%f", n_f);
x2_f = fscanf(fin, "%f", n_f);
y2_f = fscanf(fin, "%f", n_f);
fclose(fin);

fin = fopen('sss.out')
n_s = fscanf(fin, "%i", 1);
x1_s = fscanf(fin, "%f", n_s);
y1_s = fscanf(fin, "%f", n_s);
x2_s = fscanf(fin, "%f", n_s);
y2_s = fscanf(fin, "%f", n_s);
fclose(fin);

fin = fopen('ttt.out')
n_t = fscanf(fin, "%i", 1);
x1_t = fscanf(fin, "%f", n_t);
y1_t = fscanf(fin, "%f", n_t);
x2_t = fscanf(fin, "%f", n_t);
y2_t = fscanf(fin, "%f", n_t);
fclose(fin);

figure;

plot(x1_f, y1_f, 'k', x1_s, y1_s, 'm', x1_t, y1_t, 'b', x1_f, f1(x1_f), 'g');
title('e^x, user');
grid on;

legend('3 nodes', '5 nodes', '9 nodes', 'original');

figure;

plot(x2_f, y2_f, 'k', x2_s, y2_s, 'm', x2_t, y2_t, 'b', x2_f, f2(x2_f), 'g');
title('e^x/(1 - x^2), user');
grid on;

legend('3 nodes', '5 nodes', '9 nodes', 'original');