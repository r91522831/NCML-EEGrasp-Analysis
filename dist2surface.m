function dist = dist2surface(Mcom, dCOPz, dCOPy, dFy, Fgrip)
b = Mcom;
a = dCOPz;
x1 = dCOPy;
y1 = dFy;
z1 = Fgrip;

syms x y z lamda
d = (x - x1)^2 + (y - y1)^2 + (z - z1)^2;
g = x * z - a * y - b;

e1 = diff(d, x) == lamda * diff(g, x);
e2 = diff(d, y) == lamda * diff(g, y);
e3 = diff(d, z) == lamda * diff(g, z);
e4 = g == 0;

[k1, k2, k3, ~] = solve(e1, e2, e3, e4, x, y, z, lamda);

for i = 1:length(k1)
    D(i, 1) = subs(d, [x, y, z], [vpa(k1(i)), vpa(k2(i)), vpa(k3(i))]);
end

dist = double(min(sqrt(D(arrayfun(@isreal, D)))));
end