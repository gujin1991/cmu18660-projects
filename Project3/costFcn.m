function [F, G, H] = costFcn(Z, X, y, lambda, t)
% Compute the cost function F(Z)
%
% INPUTS: 
%   Z: Parameter values
% OUTPUTS
%   F: Function value
%   G: Gradient value
%   H: Hessian value
%
% @ 2011 Kiho Kwak -- kkwak@andrew.cmu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To improve the excution speed, please program your code with matrix
% format. It is 30 times faster than the code using the for-loop.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X = M x N N is smaple number
% y = 1 x N
% W = M x 1
% kexi = N x 1
% C = 1 x 1


N = size(X, 2); 
M = size(X, 1);
W = Z(1:M)';
C = Z(M+1);
kexi = Z(M+2:N+M+1)';


term = y.*(W'*X) + C*y + kexi' - 1; %1 x N
F = sum(kexi) + lambda * (W' * W) - sum(log(term)) / t - sum(log(kexi)) / t;


grad_W = 2 * lambda * W - X * (y ./ term)' / t;
grad_C = (-1 / t) * sum(y ./ term);
grad_kexi = 1 - (1 ./ term') / t - 1 ./ (t * kexi);

G = [grad_W' grad_C grad_kexi'];

sY = y .^ 2;
sTerm = term .^ 2;


HW = zeros(M, M);
for j = 1 : M
    for k = 1 : M
        HW(j, k) = (1/t)*(X(j,:) .* X(k,:)) * (sY ./ sTerm)';
        if(j == k)
            HW(j, k) = HW(j, k) + 2 * lambda;
        end        
    end
end

%HWC = (1/t)* X * (sY ./ sTerm)';
HWC = zeros(M, 1);
for j = 1 : M
    HWC(j, 1) = (1/t)*X(j, :) * (sY ./ sTerm)';
end


%HWkexi = bsxfun(@times, X, y./sTerm) / t;
HWkexi = zeros(M, N);
for i = 1 : M
    for j = 1 : N
        HWkexi(i, j) = X(i, j) * y(j) / (t * term(j).^2);
    end
end

HC = sum(sY ./ sTerm) / t;
% HC1 = (sY * (sTerm .^ (-1))') / t;


%Hkexi = diag(1./sTerm + 1./(t .* (kexi'.^2)));
Hkexi = zeros(N, N);
for i = 1 : N
    for j = 1 : N
        if(i == j)
            Hkexi(i, j) = (1/(term(i)^2) + 1/(kexi(i)^2))/t;
        end
    end
end

HCkexi = y ./ sTerm * (1 / t); %1 x N


H = [HW HWC HWkexi; HWC' HC HCkexi; HWkexi' HCkexi' Hkexi];
end



     





