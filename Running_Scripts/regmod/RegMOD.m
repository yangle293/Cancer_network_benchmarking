function f = RegMOD(d, obs, dfm, C)

% Calculate active score using RegMOD method.
%
% Input:
% d : the adjacency matrix of network. If node i and node j are linked,
%      d(i,j) = 1 for unweighted case and d(i,j) is equal to a weight for
%      weighted case.
% obs : observed active score of each node which is calculated by SNR, 
%         t-statistic or others.
%
% Output:
% f : predicted active score of each node.

n = length(d);
% deg = sum(d);
% Lapd = -1*diag(deg) + d;
% dfm = padm(tau*Lapd);    

trainm = [(1:n)', dfm];
trainy= obs;
model = svmtrain(trainy, trainm, ['-s 3 -t 4 -c ', num2str(C)] );    
[f, accuracy, prob_estimates] = svmpredict(obs, [(1:n)', dfm], model);

