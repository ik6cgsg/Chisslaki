digits(11);

x = -2 : 0.1 : 3;
poly = @(x) x .^ 3 - 2 * x .^ 2 - 4 * x + 7;
poly_der = @(x) 3 * x .^ 2 - 4 * x - 4;
y = poly(x);

X1 = -1.9354;
X2 = 1.4626;
X3 = 2.4728;

X1_p = vpa(fzero(@(x) poly(x), -2))
X2_p = vpa(fzero(@(x) poly(x), 1.5))
X3_p = vpa(fzero(@(x) poly(x), 2.5))

Eps = 0.01;

grid on;

line(x, y);
xlabel('x');
ylabel('y');
title('polynom');

hold on;

line(x, poly_der(x), 'color', 'k');

poly_phi_1 = @(x) x - 1 / 14.7 * poly(x);
line(x, poly_phi_1(x), 'color', 'r');

poly_phi_der_1 = @(x) 1 - 1 / 14.7 * poly_der(x);
line(x, poly_phi_der_1(x), 'color', 'm');
max_1 = fminbnd(@(x) -abs(poly_phi_der_1(x)), X1 - Eps, X1 + Eps);
plot(max_1, poly_phi_der_1(max_1), '-.k*');

poly_phi_2 = @(x) x + 1 / 3.6 * poly(x);
line(x, poly_phi_2(x), 'color', [0.3, 0.6, 0.1], 'marker', '*');

poly_phi_der_2 = @(x) 1 + 1 / 3.6 * poly_der(x);
line(x, poly_phi_der_2(x), 'color', [0.3, 0.3, 0.2], 'marker', '*');
max_2 = fminbnd(@(x) -abs(poly_phi_der_2(x)), X2 - Eps, X2 + Eps);
plot(max_2, poly_phi_der_2(max_2), '-go');

poly_phi_3 = @(x) x - 1 / 4.7 * poly(x);
line(x, poly_phi_3(x), 'color', [1, 0.3, 0.3], 'linestyle', '-.');

poly_phi_der_3 = @(x) 1 - 1 / 4.7 * poly_der(x);
line(x, poly_phi_der_3(x), 'color', [0.2, 0.6, 0.7], 'linestyle', '-.');
max_3 = fminbnd(@(x) -abs(poly_phi_der_3(x)), X3 - Eps, X3 + Eps);
plot(max_3, poly_phi_der_3(max_3), ':md');

legend('polynom', 'polynom der', 'poly phi 1', 'poly phi 1 der', 'max 1', 'poly phi 2', 'poly phi 2 der', 'max 2', 'poly phi 3', 'poly phi 3 der', 'max 3');

figure
grid on;

x = -5 : 0.1 : 5;
trans = @(x) 2 * x .* sin(x) - cos(x);
y = trans(x);

X1 = -0.6532;
X2 = 0.6532;
X3 = 3.2923;

X1_t = vpa(fzero(@(x) trans(x), -0.5))
X2_t = vpa(fzero(@(x) trans(x), 0.7))
X3_t = vpa(fzero(@(x) trans(x), 3))

Eps = 0.01;

line(x, y);
xlabel('x');
ylabel('y');
title('trans');

hold on;

trans_der = @(x) 3 .* sin(x) + 2 * x .* cos(x);
line(x, trans_der(x), 'color', 'k');

%legend('transcendent', 'transcendent derivative');

trans_phi_1 = @(x) x + 2 / 5.7 * trans(x);
line(x, trans_phi_1(x), 'color', 'r');

trans_phi_der_1 = @(x) 1 + 2 / 5.7 * trans_der(x);
line(x, trans_phi_der_1(x), 'color', 'm');
max_1 = fminbnd(@(x) -abs(trans_phi_der_1(x)), X1 - Eps, X1 + Eps);
plot(max_1, trans_phi_der_1(max_1), '-.k*');

trans_phi_2 = @(x) x - 2 / 5.7 * trans(x);
line(x, trans_phi_2(x), 'color', [0.3, 0.6, 0.1], 'marker', '*');

trans_phi_der_2 = @(x) 1 - 2 / 5.7 * trans_der(x);
line(x, trans_phi_der_2(x), 'color', [0.3, 0.3, 0.2], 'marker', '*');
max_2 = fminbnd(@(x) -abs(trans_phi_der_2(x)), X2 - Eps, X2 + Eps);
plot(max_2, trans_phi_der_2(max_2), '-go');

trans_phi_3 = @(x) x + 2 / 14.0 * trans(x);
line(x, trans_phi_3(x), 'color', [1, 0.3, 0.3], 'linestyle', '-.');

trans_phi_der_3 = @(x) 1 + 2 / 14.0 * trans_der(x);
line(x, trans_phi_der_3(x), 'color', [0.2, 0.6, 0.7], 'linestyle', '-.');
max_3 = fminbnd(@(x) -abs(trans_phi_der_3(x)), X3 - Eps, X3 + Eps);
plot(max_3, trans_phi_der_3(max_3), ':md');

legend('trans', 'trans der', 'trans phi 1', 'trans phi 1 der', 'max 1', 'trans phi 2', 'trans phi 2 der', 'max 2', 'trans phi 3', 'trans phi 3 der', 'max 3');

figure
grid on;

x = 1.3 : 0.01 : 2;

XP = 1.4626;
XP0 = 2;
XL = 1;
XR = 2;

title('Simple iterations for polynom');

line(x, x, 'color', 'b');
line(x, poly_phi_2(x), 'color', 'k');
xk = XP0;
arrx = ones(8, 1);
arry = arrx;

for i = 1 : 4
    xk_1 = xk;
    xk = poly_phi_2(xk);
    arrx(2 * i - 1) = xk_1;
    arrx(2 * i) = xk_1;
    arry(2 * i - 1) = xk_1;
    arry(2 * i) = xk;
end
line(arrx, arry, 'color', 'r');

legend('y = x', 'y = phi(x)',...
       '4 simple iterations steps', 2);
   
%{   
for i = 1 : 4
   tmp = (XR + XL) / 2;
   xh = tmp;
   yh1 = poly(tmp) - 1;
   yh2 = poly(tmp) + 1;
   line([xh xh], [yh1 yh2], 'color', 'g');
   if (poly(tmp) * poly(XR) <= 0)
       XL = tmp;
   else if (poly(XL) * poly(tmp) <= 0)
       XR = tmp;
       end
   end
end
%}

figure
grid on;

%x = 0 : 0.01 : 1.7;
x = 0 : 0.01 : 3;

title('Simple iterations for transcendent');

XT = 0.6532;
XT0 = 2.6;
XL = 0;
XR = 1.5;

%line(x, trans(x), 'linestyle', ':');
line(x, x, 'color', 'b');
line(x, trans_phi_2(x), 'color', 'k');

xk = XT0;
%xk = trans_phi_2(xk);

for i = 1 : 4
    xk_1 = xk;
    xk = trans_phi_2(xk);
    arrx(2 * i - 1) = xk_1;
    arrx(2 * i) = xk_1;
    arry(2 * i - 1) = xk_1;
    arry(2 * i) = xk;%trans_phi_2(xk);
end
line(arrx, arry, 'color', 'r');

legend('y = x', 'y = phi(x)',...
       '4 simple iterations steps', 0);

%{
for i = 1 : 4
   xh = tmp;
   yh1 = trans(tmp) - 1;
   yh2 = trans(tmp) + 1;
   line([xh xh], [yh1 yh2], 'color', 'g');
   tmp = (XR + XL) / 2;
   if (trans(tmp) * trans(XR) <= 0)
       XL = tmp;
   else if (trans(XL) * trans(tmp) <= 0)
       XR = tmp;
       end
   end
end
%}

%legend('trans', 'trans phi', '4 simple iterations steps', '4 half divisoins iterations');
