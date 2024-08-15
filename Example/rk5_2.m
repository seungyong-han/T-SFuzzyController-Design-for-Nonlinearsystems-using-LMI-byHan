function dx=rk5_2(x,u,T)
k1=Inverted_pendulum2(x,u)*T;
k2=Inverted_pendulum2(x+k1*0.5,u)*T;
k3=Inverted_pendulum2(x+k2*0.5,u)*T;
k4=Inverted_pendulum2(x+k3,u)*T;
dx=x + ((k1+k4)/6+(k2+k3)/3);