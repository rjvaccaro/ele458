function dsm(phi,gamma,K)
%DSM	Digital Stability margins.
%	dsm(phi,gamma,K) computes and prints the upper and lower gain margins,
%	as well as the phase margin(s) for the loop transfer function
%	described by the discrete-time state-space model (phi,gamma,K).
%	The maximum margins the function searches for are -30db lower gain 
%	margin, 30db upper gain margin, and 120 degree phase margin.

%  R.J. Vaccaro 5/94

ugm(phi,gamma,K);
fprintf('\n')
lgm(phi,gamma,K);
fprintf('\n')
phm(phi,gamma,K);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ugm(A,B,L);
%UGM	Upper gain  margin.
%	ugm(A,B,L) computes and prints the upper gain margin(s) for the 
%	loop transfer function described by the state space matrices A,B,L. 
%	The maximum margin the function searches for is 30 dB.

%  R.J. Vaccaro 5/94

UGmax=31;
[n,p]=size(B);
for loop = 1:p
  bv=B(:,loop);
  t1=1;
  t2=2;
  skip=0;
  B(:,loop)=t2*bv;
  while skip==0 & max(abs(eig(A-B*L)))<1
    if t2>UGmax
	skip=1;
	ug=t2;
    end
    t2=2*t2;
    B(:,loop)=t2*bv;
  end
  if skip==0
    while 20*log10(t2/t1)>0.001
    	tt=(t1+t2)/2;
	B(:,loop)=tt*bv;
	if  max(abs(eig(A-B*L)))>1
	  t2=tt;
	else
	  t1=tt;
	end
    end
    ug=tt;
  end
  ugm(loop)=20*log10(ug);
  fprintf('Upper gain margin for input #%g is %g dB\n',loop,...
					round(100*ugm(loop))/100);
  B(:,loop)=bv;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lgm(A,B,L);
%LGM    Lower gain  margin.
%       lgm(A,B,L) computes and prints the lower gain margin(s) for the
%       loop transfer function described by the state space matrices A,B,L.
%       The maximum margin the function searches for is -30 dB.

%  R.J. Vaccaro 5/94

LGmin=1/31;
[n,p]=size(B);
for loop = 1:p
  bv=B(:,loop);
  t1=1;
  t2=1/2;
  skip=0;
  B(:,loop)=t2*bv;
  while skip==0 & max(abs(eig(A-B*L)))<1
    if t2<LGmin
	skip=1;
	lg=t2;
    end
    t2=t2/2;
    B(:,loop)=t2*bv;
  end
  if skip==0
    while 20*log10(t1/t2)>0.001
    	tt=(t1+t2)/2;
	B(:,loop)=tt*bv;
	if  max(abs(eig(A-B*L)))>1
	  t2=tt;
	else
	  t1=tt;
	end
    end
    lg=tt;
  end
  lgm(loop)=20*log10(lg);
  fprintf('Lower gain margin for input #%g is %g dB\n',loop,...
					round(100*lgm(loop))/100);
  B(:,loop)=bv;
end
%
%  END OF LGM.M 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phm(A,B,L);
%PHM    Phase margin.
%       phm(A,B,L) computes and prints the phase margin(s) for the
%       loop transfer function described by the state space matrices A,B,L.
%       The maximum margin the function searches for is 120 degrees.

%  R.J. Vaccaro 5/94

j=sqrt(-1);
PMmax=120;
[n,p]=size(B);
cf=pi/180;
for loop = 1:p
  bv=B(:,loop);
  t1=0;
  t2=16;
  skip=0;
  B(:,loop)=exp(-j*t2*cf)*bv;
  while skip==0 & max(abs(eig(A-B*L)))<1
    if t2>PMmax
	skip=1;
	pm=t2;
    end
    t2=2*t2;
    B(:,loop)=exp(-j*t2*cf)*bv;
  end
  if skip==0
    while t2-t1 > 0.5
    	tt=(t1+t2)/2;
	B(:,loop)=exp(-j*tt*cf)*bv;
	if  max(abs(eig(A-B*L)))>1
	  t2=tt;
	else
	  t1=tt;
	end
    end
    pm=tt;
  end
  phm(loop)=pm;
  fprintf('Phase margin for input #%g is %g degrees\n',loop,...
					round(phm(loop)));
  B(:,loop)=bv;
end
%
%  END OF PHM.M 
