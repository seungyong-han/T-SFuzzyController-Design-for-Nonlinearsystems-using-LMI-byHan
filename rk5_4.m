function dx=rk5_4(x,u,T)
k1=ND(x,u)*T;
k2=ND(x+k1*0.5,u)*T;
k3=ND(x+k2*0.5,u)*T;
k4=ND(x+k3,u)*T;
dx=x + ((k1+k4)/6+(k2+k3)/3);