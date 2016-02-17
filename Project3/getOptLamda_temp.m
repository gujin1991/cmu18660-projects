function [optSolution, optLamda] = getOptLamda_temp(X, Y, setPara)
% Get the optimal lamda
%
% INPUTS:
%   X(MxN) : trData(i,j) is the i-th feature from the j-th trial
%   Y(Nx1): trData(j) is the label of the j-th trial (1 or -1)
%   setPara : Initialized parameters
%            setPara.t      
%            setPara.beta   
%            setPara.Tmax   
%            setPara.tol    
%            setPara.W      
%            setPara.C      
%
% OUTPUTS:
%   optiLamda: Optimal lamda value 
%
% @ 2011 Kiho Kwak -- kkwak@andrew.cmu.edu
fea_num = size(X, 1);
W = setPara.W;
C = setPara.C;
tol = setPara.tol;
Tmax = setPara.Tmax;
optLamda = 0;
beta = setPara.beta;
error = Inf;
fold_num = 5;
lambda_list = [0.01 1 100 10000];

c1 = X(:,1:100);
c2 = X(:,101:200);
label1 = Y(1:100);
label2 = Y(101:200);

for l = 1 : size(lambda_list, 2)
    lambda = lambda_list(l);
    localError = 0;
    dis = size(c1, 2) / fold_num;
    for i = 1 : fold_num
        testInx = (i -1)*dis+1 : i*dis;
        testData = [c1(:,testInx) c2(:,testInx)];
        trainInx = setdiff(1:100, testInx);
        trainData = [c1(:,trainInx) c2(:,trainInx)];
        testLab = [label1(:,testInx) label2(:,testInx)];
        trainLab = [label1(:,trainInx) label2(:,trainInx)];
        
        zeta = zeros(size(trainData,2), 1);
        for j = 1 : size(trainData,2)
            zeta(j,1) = max(1 - trainLab(j)*(W' * trainData(:, j)+C), 0) + 0.001;
        end
        
        init_Z = [W', C, zeta'];
        %test_value = trainLab.*(W'*trainData) + C + zeta;
        t = setPara.t;
        while (t <= Tmax)
            [localSol, err] = solveOptProb_NM_temp(@costFcn, init_Z, tol, trainData, trainLab, lambda, t);
            init_Z = localSol;
            t = t * beta;
        end
        localW = localSol(1:fea_num);
        localC = localSol(fea_num+1);
        testVal = (localW * testData + localC);
        localError = localError + sum(((testVal .* testLab) < 0));
    end
    
    if(localError < error)
        error = localError;
        optLamda = lambda;
        optSolution = localSol;
    end
end
