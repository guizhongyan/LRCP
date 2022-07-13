function alpha = Updata_w(S,A)
%the code is written by Kun zhan in Graph learning for multiview clustering(MVGL).
[n,~,nv] = size(A);
% mn = m*n;
% % s = reshape(S, [mn 1]);
alpha = zeros(nv,n);
for k = 1:n
    s = S(:,k);
    Z = zeros(n,nv);
    for i = 1:nv
        a = A(:,k,i);
    %     a = reshape(A(:,:,i), [mn 1]);
        Z(:,i) = s - a;
    end
    C = Z'*Z;
    C = C + eye(nv)*eps*trace(C);
    W = C\ones(nv,1);
    alpha(:,k) = W/sum(W);
end