function [ Temperature ] = thermalsimGauss( p, mediumX, mediumY, leftBound, rightBound, topBound, bottomBound )
%THERMALSIMGAUSS solves the 2D steady state thermal problem using Gaussian
%elimination
%   INPUT:
%   p:  discretized power density
%   mediumX:    x-dimension of the medium
%   mediumY:    y-dimension of the medium
%	leftBound:	Temperature at the left boundary (x=0), leftBound(j) means
%	the temperature at T(0,j)
%	rightBound:	Temperature at the right boundary (x=N+1)
%	topBound:	Temperature at the top boundary (y=M+1)
%	bottomBound:	Temperature at the bottom boundary (y=0)
%
%   OUTPUT:
%   Temperature: solved thermal map
tic;
N = size(topBound, 1);
M = size(leftBound, 1);
deltaX = mediumX / N;
deltaY = mediumY / M;
sX = deltaX ^ 2;
sY = deltaY ^ 2;
sumXY = sX + sY;
mulXY = sX * sY;
four = (2 * sumXY / mulXY) * eye(M, M);
one = (-1 / sX) * eye(M, M);
zero = zeros(M, M);
k = 157;


for i = 1 : M
    for j = 1 : M
        if abs(i - j) == 1
            four(i, j) = -1 / (sY);
        end
    end
end

p = p ./ k;
for i = 1 : N
    for j = 1 : M
        if(i == 1)
            p(i, j) = p(i, j) + leftBound(j,:) / sY;
        end
        if(i == N)
            p(i, j) = p(i, j) + rightBound(j,:) / sY;
        end
        if(j == 1)
            p(i, j) = p(i, j) + bottomBound(i,:) / sX;
        end
         if(j == M)
            p(i, j) = p(i, j) + topBound(i,:) / sX;
         end
    end
end

b = p(1, :)';
for i = 2 : N
    b = [b; p(i, :)'];
end


for i = 1 : N
    AA = zeros(M, M);
    for j = 1 : N
        if j == 1
            if i == 1
            AA = four;
            else if i == 2
                  AA = one;
                else
                    AA = zero;
                end
            end           
        else if i == j
                AA = [AA, four];
            else if abs(i - j) == 1
                    AA = [AA, one];
                else
                    AA = [AA, zero];
                end
            end
        end
    end
    if i == 1
        A = AA;
    else
        A = [A; AA];
    end
end

trueT = A \ b;
B = [A b]; 
n = length(b);
RA = rank(A);
RB = rank(B);
if abs(RB - RA) > 0
    disp('There is no answer we can find. Because the rank of coeffcient matrix A is not equal to the rank of augmented matrix B.\n')
    return
end

if RA == RB
    if RA == n 
        X = zeros(n,1);
        for i = 1 : n - 1  %change the matrix to be a upper triangular matrix                                    
              for j = i + 1 : n 
                    if B(j, i) ~= 0
                    coe = B(j,i) / B(i,i);   
                    B(j,i:n+1) = B(j,i:n+1) - coe * B(i,i:n+1);
                    end
              end
        end
        
        % solve the equatuon set from Xn to X1
        b = B(1:n, n+1);
        A = B(1:n, 1:n);
        X(n) = b(n) / A(n,n);
        for i = n - 1 :-1: 1
            X(i) = (b(i) - A(i, i+1:n ) * X(i+1:n)) / A(i,i);
        end
        
       
    else
        disp('There are more than one solutions because rank of coefficient matrix A is smaller than N.\n')
        return
    end
end



error = sqrt(sum(trueT - X) .^2 / sum(trueT .^2))
Temperature = reshape(X, M, N);
Temperature = Temperature';
toc;   


end
