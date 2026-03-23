function [delta1,delta2]=rb_regob(A,B,C,K,L,T)
%RB_REGOB       Computes input-multiplicative stability robustness bound
%               for observer-based regulator
% [delta1,delta2]=rb_regob(A,B,C,K,L,T)
%INPUTS
%  A,B,C        State-space plant model (continuous- or discrete-time)
%  K            State feedback gain matrix
%  L            Observer gain matrix
%  T            Sampling interval (Use T=0 for continuous-time system)
%
%OUTPUT
%  delta1       Unstructured input-multiplicative infinity-norm 
%               perturbation bound 
%  delta2       Unstructured input-multiplicative feedback 
%               infinity-norm perturbation bound
%
% R.J. Vaccaro 3/2016
[n,p]=size(B);
D1=zeros(p,p);
    sys1=ss([A -B*K;L*C A-L*C-B*K],[B; 0*B],[0*K -K],D1,T);
    delta1=1/norm(sys1,inf);
    sys2=ss([A -B*K;L*C A-L*C-B*K],[B; 0*B],[0*K -K],eye(p),T);
    delta2=1/norm(sys2,inf);
end
