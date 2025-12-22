function [K1,K2,delta1,delta2]=tsd(A,B,C,Aa,Ba,poles,T,algorithm);
%TSD	Tracking system design (analog or digital) with specified algorithm
%	[K1,K2,delta1,delta2]=tsd(A,B,C,Aa,Ba,poles,T,algorithm)
%
% INPUTS
%   (A,B,C)         Plant model (continuous or discrete time)
%   (Aa,Ba)         Additional dynamics model
%   poles           Desired closed-loop pole locations for the
%                   design model (cascade of plant and additional
%                   dynamics)
%   T               Sampling interval (use 0 for continuous-time systems)
%   algorithm       'fbg' or 'place' or 'rfbg'
%
% OUTPUTS
%
%   K1 and K2       Feedback gains, K1 for plant state variables, K2 for
%                   additional dynamics
%   delta1, delta2  Stability robustness bounds
%
%  R.J. Vaccaro  3/2016
[n,p]=size(B);
[q,m]=size(Ba);
Ad=[A zeros(n,q);Ba*C Aa];
Bd=[B;zeros(q,p)];
if strcmp(algorithm ,'rfbg')
    eval(strcat('[Kd,delta1,delta2]=',algorithm,'(Ad,Bd,poles,T);'))
    K1=Kd(:,1:n);
    K2=Kd(:,n+1:n+q);  
elseif strcmp(algorithm ,'place') | strcmp(algorithm ,'fbg')
    eval(strcat('Kd=',algorithm,'(Ad,Bd,poles);'))
    K1=Kd(:,1:n);
    K2=Kd(:,n+1:n+q);
    [n,p]=size(B);
D1=zeros(p,p);
    sys=ss(Ad-Bd*Kd,Bd,Kd,D1,T);
    delta1=1/norm(sys,inf);
    sys=ss(Ad-Bd*Kd,Bd,-Kd,eye(p),T);
    delta2=1/norm(sys,inf);
end
return
