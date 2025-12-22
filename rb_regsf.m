function [delta1,delta2]=rb_regsf(A,B,K,T)
%RB_REGSF       [delta1,delta2]=rb_regsf(A,B,K,T)  Computes input-multiplicative 
%               and input-feedback stability robustness bounds for an analog or 
%		digital state-feedback regulator
%
%INPUTS
%  A,B          State-space plant model (continuous- or discrete-time)
%  K            State feedback gain matrix
%  T            Sampling interval (Use T=0 for continuous-time system)
%
%OUTPUT
%  delta1       Input-multiplicative stability robustness bound
%  delta2       Input-feedback stability robustness bound
%
%               The following guaranteed classical stability margins are 
%		available simultaneously on all plant inputs:
%
%               Upper gain margin >= 20 log10(max(1+delta1,1/(1-delta2))) dB
%               Lower gain margin <= 20 log10(min(1-delta,1/(1+delta2))) dB
%               Phase Margin >= 2*180/pi*asin(max(delta1,delta2)/2) degrees
%
% R.J. Vaccaro 11/2016
[n,p]=size(B);
D1=zeros(p,p);
    sys=ss(A-B*K,B,K,D1,T);
    delta1=1/norm(sys,inf);
    sys=ss(A-B*K,B,-K,eye(p),T);
    delta2=1/norm(sys,inf);
end
