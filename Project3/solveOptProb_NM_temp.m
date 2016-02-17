function [optSolution, err] = solveOptProb_NM_temp(costFcn,init_Z,tol, X, y, lambda, t)
% Compute the optimal solution using Newton method
%
% INPUTS:
%   costFcn: Function handle of F(Z)
%   init_Z: Initial value of Z
%   tol: Tolerance
%
% OUTPUTS:
%   optSolution: Optimal soultion
%   err: Errorr
%
% @ 2011 Kiho Kwak -- kkwak@andrew.cmu.edu


fea_num = size(X,1);
data_num = size(X,2);

Z = init_Z;  
[F, G, H] = feval(costFcn,Z,X,y,lambda,t);
deltaZ = -inv(H) * G';
err = -(G * deltaZ);


while (err/2) > tol    
    s = 1;
    while(true)  
        temp = Z + s*deltaZ'; 
        localW = temp(1:fea_num);
        localC = temp(fea_num+1);
        localkesi = temp(fea_num + 2:fea_num+1 + data_num);
        if(isempty((find(((localW * X) .* y + localC * y + localkesi - 1)<=0,1))) && isempty((find(localkesi <=0,1)))),
            break;
        else
            s = 0.5 * s;            
        end           
    end
    Z = Z + s*deltaZ'; 
    
    [F, G, H] = feval(costFcn,Z,X,y,lambda,t);
    deltaZ = -inv(H) * G';
    err= -(G * deltaZ);
 end
 optSolution = Z;
end






