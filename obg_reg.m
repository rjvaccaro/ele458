function [L,delta1,delta2]=obg_reg(A,B,C,K,poles,T)
%OBG_REG  Observer gain calculation for observer-based regulator
%[L,delta1,delta2] =obg_reg(A,B,C,K,poles,T)
%
%INPUTS
%  A,B,C        State-space plant model (continuous- or discrete-time)
%  K            State-feedback gain matrix
%  poles        Vector of desired observer poles
%  T            Sampling interval (Use T=0 for continuous-time system)
%
%OUTPUT
%  L            Observer gain matrix. Note: eig(A-L*C)=poles
%  delta1       Unstructured input-multiplicative infinity-norm 
%               perturbation bound 
%  delta2       Unstructured input-multiplicative feedback 
%               infinity-norm perturbation bound
%
% R.J. Vaccaro 3/2016
%
% Functions all_pos_combs.m and fmin_V2 were written by Edgar Ponce 
%          as part of an ELE 598 project in the spring of 2024.
%
[m,n]=size(C);
[W1,W2,cvo] = fmin_V2(A,C,poles);
nt=n*(m-1); % number of theta parameters needed

theta=pi/4*ones(nt,1);
tfun=[.01 .01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01];
tx=[0.01 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001];
fprintf('\ndelta1   delta2 \n')
for LOOP=1:5
    options=optimset('display','off','TolX',tx(LOOP),'TolFun',tfun(LOOP),...
        'MaxFunEvals',5000);
[theta,f]=fminsearch(@cost,theta,options,A,B,C,K,cvo,W1,W2,T);
L=formL(W1,W2,cvo,theta,n,m)';
[delta1,delta2]=rb_regob(A,B,C,K,L,T);
fprintf('%1.4f   %1.4f\n',delta1,delta2)
end
L=formL(W1,W2,cvo,theta,n,m)';
%--------------------------------------------------------------------------
function [L,cond_flg]=formL(V1,V2,cv,theta,n,p)
cond_flg=0;
ncp=sum(cv); % number of complex poles/2
nrp=n-2*ncp; % number of real poles
M=[];
if nrp>0
    for k=1:nrp
        theta1=theta((k-1)*(p-1)+1:k*(p-1));
        alpha1=hyperspherical(theta1);
        M=[M V1(:,(k-1)*p+1:k*p)*alpha1];
    end
end
if ncp>0
    offset1=nrp*(p-1);
    offset2=offset1+ncp*(p-1);
    for k=1:ncp
        theta1=theta(offset1+(k-1)*(p-1)+1:offset1+k*(p-1));
        alpha1=hyperspherical(theta1);
        theta2=theta(offset2+(k-1)*(p-1)+1:offset2+k*(p-1));
        theta2=[theta2;-sum(theta2)]; % RJV
        alpha2=alpha1.*exp(1i*theta2);
        X=V2(:,(k-1)*p+1:k*p)*alpha2;
        M=[M real(X) imag(X)];
    end
end
if cond(M(1:n,:))>1.e9
    L=[];
    cond_flg=1;
else
    L=M(n+1:end,:)/M(1:n,:);
end
%--------------------------------------------------------------------------
function h=hyperspherical(theta)
%theta is a vector containing n-1 angles (radians)
%h is a vector containing n Cartesian coordinates of
%a point on a unit hypersphere
p=length(theta);
h=[cos(theta);1];
for k=1:p
    h=h.*[ones(k,1);sin(theta(k))*ones(p-k+1,1)];
end
%--------------------------------------------------------------------------
function f=cost(theta,A,B,C,K,cvo,V1,V2,T)
[m,n]=size(C);
[L,cond_flg]=formL(V1,V2,cvo,theta,n,m);
L=L';
if ~cond_flg
[n,p]=size(B);
D1=zeros(p,p);
    sys1=ss([A -B*K;L*C A-L*C-B*K],[B; 0*B],[0*K -K],D1,T);
    f1=norm(sys1,inf);
    delta1 = 1/f1;
    sys2=ss([A -B*K;L*C A-L*C-B*K],[B; 0*B],[0*K -K],eye(p),T);
    f2=norm(sys2,inf);
    delta2 = 1/f2;
%     f=f1+100*f2;
    f = 1/(delta1+delta2);
else
    f=1.e10;
end
%--------------------------------------------------------------------------
function pole_combs = all_pos_combs(poles)
n = 2^sum((imag(poles) ~= 0)/2);
pole_combs(1, :) = poles;
cp_indices = find(imag(poles) ~= 0); 
combs = dec2bin(0:n-1)-'0';

for i = 2:n
    ra_poles = poles;
    for j = 1:length(cp_indices)/2 
        if combs(i,j) == 1
            ra_poles([cp_indices(2*j-1), cp_indices(2*j)]) = ...
                poles([cp_indices(2*j), cp_indices(2*j-1)]);
        end
    end
    pole_combs(i, :) = ra_poles;
end
%--------------------------------------------------------------------------
function [V1,min_V2,cvo] = fmin_V2(A,C,poles)
pole_combs = all_pos_combs(poles);
min_condition_number = Inf;
for i = 1:2^sum((imag(poles) ~= 0)/2)
    poles = pole_combs(i, :);
    cvo = ~imag(poles) == 0;
    ind = find(cvo);

    if ~isempty(ind)
        ind = ind(1:2:end);
        cvo(ind) = [];
        poles(ind) = [];
    end

    [V1,V2]=obasis2(A',C',cvo,poles);
    condition_V2 = cond(V2);
        
    if condition_V2 < min_condition_number
        min_condition_number = condition_V2;
        min_V2 = V2;
    end
end
%--------------------------------------------------------------------------
function [V1,V2]=obasis2(A,B,cv,poles)
[n,p]=size(B);
ncv=length(cv);
scv=sum(cv);
V1=[];
V2=[];
if n-2*scv>0
    V1=zeros(n+p,(n-2*scv)*p); % basis vectors for real-valued poles
end
if scv>0
    V2=zeros(n+p,scv*p); % basis vectors for complex-valued poles
end
c1=0;
c2=0;
for k=1:length(cv)
    X=null([(poles(k)*eye(n)-A) B]);
    if ~cv(k)
        c2=c2+1;
        V1(:,(c2-1)*p+1:c2*p)=X;
    else
        c1=c1+1;
        V2(:,p*(c1-1)+1:p*c1)=orth(X); %RJV
    end
end
%--------------------------------------------------------------------------
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

