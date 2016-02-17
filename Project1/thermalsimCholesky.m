function [ Temperature ] = thermalsimCholesky( p, mediumX, mediumY, leftBound, rightBound, topBound, bottomBound )
%THERMALSIMCHOLESKY solves the 2D steady state thermal problem using
%Cholesky factorization
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

p = p / k;
for i = 1 : N
    for j = 1 : M
        if(i == 1)
            p(i, j) = p(i, j) + leftBound(j,:) / sX;
        end
        if(i == N)
            p(i, j) = p(i, j) + rightBound(j,:) / sX;
        end
        if(j == 1)
            p(i, j) = p(i, j) + bottomBound(i,:) / sY;
        end
         if(j == M)
            p(i, j) = p(i, j) + topBound(i,:) / sY;
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
n = size(A, 1);
L = zeros(size(A));
for i = 1 : n
    for j = 1 : n
        
        if(i == j)
            L(i, j) = sqrt(A(i, j) - (L(j, 1:j -1) * L(j, 1:j -1)'));
        end
        if(i > j) 
            L(i, j) = (A(i, j) - (L(i, 1:j -1) * L(j, 1:j -1)')) / L(j, j);
        end
    end
end

y = zeros(n, 1);
X = zeros(n, 1);
y(1) = b(1) / L(1, 1);
for i = 2 : n
   y(i) = (b(i) - L(i, 1:i-1 ) * y(1:i-1)) / L(i,i);
end

L = L';
X(n) = y(n) / L(n,n);
for i = n - 1 :-1: 1
    X(i) = (y(i) - L(i, i+1:n ) * X(i+1:n)) / L(i,i);
end

error = sqrt(sum(trueT - X) .^2 / sum(trueT .^2))
Temperature = reshape(X, M, N);
Temperature = Temperature';
toc;

end
