function [result C]=clustering(S, cls_num, gt)

[C] = SpectralClustering(S,cls_num);
[result] = Clustering8Measure(gt,C);      

end
