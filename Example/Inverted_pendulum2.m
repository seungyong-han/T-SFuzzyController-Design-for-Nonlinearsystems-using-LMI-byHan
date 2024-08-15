function dx=Inverted_pendulum2(x,u)
global f0 f1 m M l J g
dx(1,1) = x(2);
a1 = -f1*(M+m)*x(2) - (m^2*l^2*x(2)^2)*sin(x(1))*cos(x(1))+...
    f0*m*l*x(4)*cos(x(1)) + (M+m)*m*g*l*sin(x(1)) - m*l*cos(x(1))*u;
a2 = (M+m)*(J+m*l^2) - (m^2*l^2*cos(x(1))^2);

dx(2,1) = a1/a2;
dx(3,1) = x(4);

b1 = f1*m*l*x(2)*cos(x(1)) + (J+m*l^2)*m*l*x(2)^2*sin(x(1))+...
    -f0*(J+m*l^2)*x(4)-(m^2*g*l^2)*sin(x(1))*cos(x(1)) + (J+m*l^2)*u;
b2 = (M+m)*(J+m*l^2) - (m^2*l^2*cos(x(1))^2);
dx(4,1) = b1/b2;

