function z = Function3(xv)

x = xv(1);
y = xv(2);
z = 100*(y - x.^2).^2 + (1-x).^2 + 2;