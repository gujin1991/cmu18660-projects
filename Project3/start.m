load('feaSubEImg.mat');
%load('feaSubEOvert.mat');
fold_num = 6;
c1 = [class{1}];
c2 = [class{2}];
data_num = size(c1, 2);
fea_num = size(c1, 1);
label1 = ones(1, data_num);
label2 = -ones(1, data_num);

res = cell(fold_num, 1);
dis = size(c1, 2) / fold_num;
Accuracy = zeros(fold_num, 1);
setPara.W = ones(fea_num, 1);
setPara.C = 0;
setPara.t = 1;
setPara.tol = 0.000001;
setPara.Tmax = 1000000;
setPara.beta = 15;
for i = 1 : fold_num
    tic;
    testInx = (i -1)*dis+1 : i*dis;
    testData = [c1(:,testInx) c2(:,testInx)];
    trainInx = setdiff(1:size(c1, 2), testInx);
    trainData = [c1(:,trainInx) c2(:,trainInx)];
    testLab = [label1(:,testInx) label2(:,testInx)];
    trainLab = [label1(:,trainInx) label2(:,trainInx)];
    [optSolution, optLambda] = getOptLamda_temp(trainData, trainLab, setPara);    
    res{i} = optSolution;
    W = optSolution(1:fea_num);
    C = optSolution(fea_num + 1);  
    val = (W * testData + C);
    Accuracy(i) = sum((val .* testLab) > 0) / size(testData, 2);
    disp(Accuracy(i));
    show_chanWeights(abs(W));
    toc;
 
end
