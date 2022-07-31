function [S] = Updata_A(A0,c,islocal,H,r1)

    if nargin < 3
        islocal = 1;
    end;

    A0 = A0-diag(diag(A0));
    num = size(A0,1);
    dist=H;
    S = zeros(num);
    for i=1:num
        a0 = A0(i,:);
        if islocal == 1
            idxa0 = find(a0>0);
        else
            idxa0 = 1:num;
        end;
        ai = a0(idxa0);
        di = dist(i,idxa0);
        ad=ai+(r1/(r1+1))*di;
        S(i,idxa0) = EProjSimplex_new(ad);
    end;

