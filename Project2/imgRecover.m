function imgOut = imgRecover(imgIn, blkSize, numSample)
% Recover the input image from a small size samples
%
% INPUT:
%   imgIn: input image
%   blkSize: block size
%   numSample: how many samples in each block
%
% OUTPUT:
%   imgOut: recovered image
%
% @ 2011 Huapeng Zhou -- huapengz@andrew.cmu.edu
imgOut = zeros(size(imgIn));
M = 20;  %iteration number
lambda = 3:3:numSample + 10;
Size = blkSize * blkSize;
T = zeros(Size, Size);

%get DCT matrix T
for i = 1 : Size
    y = ceil(i / blkSize);
    x = i - (y - 1) * blkSize;
    for j = 1 : Size
        v = ceil(j / blkSize);
        u = j - (v - 1) * blkSize;
        if(u == 1)
            au = sqrt(1 / blkSize);
        else
            au = sqrt(2 / blkSize);
        end
        
        if(v == 1)
            bv = sqrt(1 / blkSize);
        else
            bv = sqrt(2 / blkSize);
        end
        
        MulXU = pi * (2 * x - 1) * (u - 1) / (2 * blkSize);
        MulYV = pi * (2 * y - 1) * (v - 1) / (2 * blkSize);
        T(i, j) = au * bv * cos(MulXU) * cos(MulYV); 
    end
end


for i = 1 : size(imgIn, 1) / blkSize
    for j = 1 : size(imgIn, 2)/ blkSize   
		%get one block to process
        block = imgIn(blkSize*(i-1) + 1:blkSize*(i-1)+blkSize, blkSize*(j-1) + 1:blkSize*(j-1)+blkSize);
        [B, idx] = datasample(reshape(block, 1, Size), numSample, 'Replace', false);  %sampling data
        B = B';
        A = T(idx, :);
        
        
        m = floor(numSample / 6);
        error = zeros(size(lambda, 2), 1);
		
		%cross validation
        for count = 1 : M       
            for lam = 1 : size(lambda, 2) 
                [trainB, idxT] = datasample(B, numSample - m, 'Replace', false);
                trainA = A(idxT, :);
                alpha = OMP(trainA, trainB, lambda(:, lam));
                C = T * alpha;
                idxTest = zeros(1, size(A, 1));
                idxTest(idxT) = 1;
                idxTest = ~idxTest;
                testData = B(idxTest,:);
                predictData = C(idx, :);
                predictData = predictData(idxTest, :);
                dif = norm(testData - predictData);
                if((count ~= 1 )&& (dif > 10 * error(lam, :) / (count - 1)))
                    error(lam, :) = error(lam, :) + 10 * error(lam, :) / (count - 1);
                else
                    error(lam, :) = error(lam, :) + dif;
                end
            end           
        end
		error = error ./ M;  %find minimum error to determine the best lambda
        [~, minidx] = min(error);
        
		%do the recover process again with best lambda
        alpha = OMP(A, B, lambda(:, minidx));
        C = T * alpha;
        recover_block = reshape(C, blkSize, blkSize);
        imgOut(blkSize*(i-1) + 1:blkSize*(i-1)+blkSize, blkSize*(j-1) + 1:blkSize*(j-1)+blkSize) = recover_block;
    end
end

%lambda(:, minidx)           %%output the best lambda if you want

end