function dsm_regob(phi,gamma,C,K,L)
%DSM_REGOB	Stability margins for digital observer-based regulator.
%	dsm_regob(phi,gamma,C,K,L) computes and prints the upper and lower gain margins,
%	as well as the phase margin(s) for the loop transfer function
%	for a digital observer-based regulator.
%	The maximum margins the function searches for are -30db lower gain 
%	margin, 30db upper gain margin, and 120 degree phase margin.
%  INPUTS:
%
%  phi,gamma,C        ZOH equivalent plant model.
%  K                  Feedback gains.
%  L                  Observer gains.
%
%  R.J. Vaccaro  11/94,11/98

[nn_,mm_]=size(phi);
phi_L = [phi zeros(nn_,mm_);L*C phi-L*C-gamma*K];
[nn_,mm_]=size(gamma);
gamma_L=[gamma;zeros(nn_,mm_)];
C_L=[zeros(mm_,nn_) K];
clear nn_ mm_
dsm(phi_L,gamma_L,C_L)
