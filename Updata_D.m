function [D] = Updata_D(G,r4,miu)
D=max(abs(G)-r4/miu,0).*sign(G);
end