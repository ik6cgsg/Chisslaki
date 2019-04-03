%%% Random
m = load('randmatr');
b = load('brand');
m1 = load('randmatroff');
b1 = load('brandoff');

digits(10);

x = m \ b;
x_rand = vpa(x)
x1 = m1 \ b;
x1_rand_off = vpa(x1)

det_rand = det(m)
cond_rand = cond(m, inf)

dx_x = vpa(norm(x1 - x, inf) / norm(x, inf))
da_a = vpa(norm(m1 - m, inf) / norm(m, inf))

x = m \ b;
vpa(x);
x1 = m \ b1;
vpa(x1);

vpa(norm(x1 - x, inf) / norm(x, inf))
vpa(norm(b1 - b, inf) / norm(b, inf))

%%% Hilbert
m = hilb(10);
b = load('bhilb');
m1 = load('hilboff');
b1 = load('bhilboff');

x = m \ b;
x_hilb = vpa(x)
x1 = m1 \ b;
x_hilb_off = vpa(x1)

det_hilb = det(m)
cond_hilb = cond(m, inf)

vpa(norm(x1 - x, inf) / norm(x, inf))
vpa(norm(m1 - m, inf) / norm(m, inf))

x = m \ b;
vpa(x)
x1 = m1 \ b;
vpa(x1)

vpa(norm(x1 - x, inf) / norm(x, inf))
vpa(norm(b1 - b, inf) / norm(b, inf))



