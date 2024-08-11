function dx=ND(x,u)
global a b
dx(1,1) = x(1)+x(2)+sin(x(3))-0.1*x(4)+(x(1)^2+1)*u;
dx(2,1) = x(1)-2*x(2);
dx(3,1) = x(1)+(x(1)^2*x(2)) -0.3*x(3);
dx(4,1) = sin(x(3)) - x(4);


