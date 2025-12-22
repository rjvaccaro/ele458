load sroots
% noise=0 gives noise-free sensors
noise=0; % set equal to 1 to simulate noisy sensors
alpha1=7.84e-3;
alpha2=9.09e-3;
beta=0.225;
g=981;
xe1=15;
ue=alpha1*sqrt(2*g*xe1)/beta;
xe2=alpha1^2/alpha2^2*xe1;
xe=[xe1;xe2];
f1=alpha1*sqrt(g/(2*xe1));
f2=alpha2*sqrt(g/(2*xe2));
A=[-f1 0;f1 -f2];
B=[beta;0];
C=[0 1];
Ts=20;
T=0.1;
[phi,gamma]=c2d(A,B,T);
spoles=s3/Ts;
zpoles=exp(T*spoles);
phia=1;
gammaa=1;
[K1,K2,delta1,delta2]=tsd(phi,gamma,C,phia,gammaa,zpoles,T,'place')
dsm_tssf(phi,gamma,C,phia,gammaa,K1,K2)
