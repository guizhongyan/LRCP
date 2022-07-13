clear all
clear memory
close all;
clc
addpath('./tools');
addpath('./data');
load('HW2sources.mat');       islocal_1 = 1;islocal_2 = 1; kk = 15; 	%final
X=data;
for i=1:length(X)
    X{i} = X{i}';
end
Y=truelabel{1};

rho = 1.2;
miu = 1.0000e-04;
max_miu = 1e6; 
maxIters =40;

tic
k=100;
truth=Y;
gt = truth;
c = max(truth);
V=length(X);
N = size(X{1},1);

Z=cell(1,V);
S=zeros(N,N);
S1=zeros(N,N,V);
SvV = S1;

for v = 1:V
    Z{v} = Updata_Sv(X{v}',c,kk, islocal_1);% COIL20.mat
    S1(:,:,v) = Z{v};
    SvV(:,:,v) = S1(:,:,v)./V;
end
S0V = sum(SvV,3);

for i=1:V
    D{i} = size(Z{i},1); % dimension of each view
end
SD = 0;
MZ = [];
for i=1:V
    SD = SD + D{i};
    MZ = [MZ;Z{i}];  %X
end



               r1 =0.1;
               r2 =10;
               r3 =1;
               r4 =0.001;
                  Q = zeros(k,SD);
                  P = zeros(k,N);

                  A = eye(k,N);
                  B = eye(k,N);
                  M = eye(N,N);
                  D = M;
                  YY = zeros(N,N);
                  S0 = S0V;

                   for it=1:maxIters

                        %------update S, W-------
                        H = P'*B;                             
                        [S]=Updata_A(S0,c,islocal_2,H,r1);                          
                        alpha = Updata_w(S,S1);
                        for i = 1:N
                            for v = 1:V        
                                Sv(:,i,v) = (1/(1+r1))*alpha(v,i).*S1(:,i,v);
                            end
                        end
                        S0 = sum(Sv,3);
                        clear H

                        %------update Q--------
                       
                       J = MZ*A';
                       [U,~,Vs] = svd(J,'econ');
                       Q = Vs*U';
                       clear J U Vs

                        %------update P--------
                        G = S*B';
                        [U,~,Vs] = svd(G,'econ');
                        P=Vs*U';
                        clear G U Vs   

                        %------update B--------
                        B = (r1 * P * S + r3* A * M')/( r1*eye(size(M)) + r3 * M * M' );

                        %------update A--------
                        A=(r2+r3)\(r2*Q*MZ+r3*B*M);

                        %------update M--------
                        M = (2 *r3 * B' * B + miu * eye(N,N))\(2* r3* B' * A + miu * D + YY);

                        %------update D--------
                       G = M - YY/miu; 
                       D = Updata_D(G,r4,miu);
                       clear G

                      % -------- Update Y-------- %           
                        YY= YY + miu*(D-M);
                        miu = min(rho*miu,max_miu);

                        Obj(it) =sum(sum((S-S0).^2)) + r1*norm((S-P'*B),'fro')^2 + r2*norm(( MZ-Q'*A),'fro')^2 + r3*norm((A-B*M),'fro')^2 + r4*norm((M),1);
                        if (it>1 && (abs(Obj(it)-Obj(it-1))/Obj(it-1)) < 10^-2)
                            break;
                        end

                   end

                  S = (S+S')/2;
                  result = clustering(S, c, gt)
  
    
                