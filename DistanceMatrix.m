function DM = DistanceMatrix(dsites,ctrs)
% Compute Distance Matrix between Evaluation/collocation points and center points 
% Inputs:
% dsites: Evaluation/collocation points
% ctrs: center points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
[N,s] = size(dsites); 
[M,s] = size(ctrs);
DM = zeros(N,M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Accumulate sum of squares of coordinate differences
for d=1:s
    DM = DM + (repmat(dsites(:,d),1,M)-repmat(ctrs(:,d)',N,1)).^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DM=DM.^.5; 
