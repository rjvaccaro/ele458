function dsm_tssf(phi,gamma,C,phia,gammaa,K1,K2)
%DSM_TSSF  Stability margins for digital state-feedback tracking system.
%	dsm_tssf(phi,gamma,C,phia,gamma,K1,K2) computes and prints the upper 
%       and lower gain margins, as well as the phase margin(s) 
%	The maximum margins the function searches for are -30db lower gain 
%	margin, 30db upper gain margin, and 120 degree phase margin.
%  INPUTS:
%
%  phi,gamma,C        ZOH equivalent plant model
%  phia,gammaa        Additional dynamics state-space model
%  K1,K2              State feedback and additional dynamics gains
%  L                  Observer gains.
%
%  R.J. Vaccaro  11/94,11/98

[n,p]=size(gamma);
[q,m]=size(gammaa);
Kd=[K1 K2];
phid=[phi zeros(n,q);gammaa*C phia];
gammad=[gamma;zeros(q,p)];
dsm(phid,gammad,Kd);

