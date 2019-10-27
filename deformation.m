% influential range ro (in mm)
ro = 16;
Ko = 2.03;
K1 = 0.438;
K2 = 0.102;

b1 = 5073;
b2 = 39.24;

ko = 909.9;
k1 = 1522;
k2 = 81.18;
figure(1)

Xo = 0.001;
X1 = 0;
X2 = 0;
a_old = 0;
a_new = 0;
v = 0;
M = 0.01;
A = 1;
B = 0.0001;
for t = 0:10000
    v = v + 0.5*(a_old + a_new);
    dXo = -sign(Xo)*abs(Xo^A*B);
    dX1 = 1.0*K1*(Xo-X1)*exp(k1*(Xo-X1))/b1;
    dX2 =  1.0*K2*(Xo-X2)*exp(k2*(Xo-X2))/b2;

    Xo = Xo + dXo;
    X1 = X1 + dX1;
    X2 = X2 + dX2;

    F = Ko*Xo*exp(ko*Xo)+K1*(Xo-X1)*exp(k1*(Xo-X1))+K2*(Xo-X2)*exp(k2*(Xo-X2));
    a_old = a_new;
    a_new = -Ko*Xo*exp(ko*Xo)/M;
    
    plot(t,F,'r',t,Xo,'b',t,X1,'g')
    hold on;
end