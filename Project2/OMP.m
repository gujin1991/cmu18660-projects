function alpha = OMP(A, B, lambda)
%step 1：initialization
F = B;
omiga = zeros(lambda, 1);
alpha = zeros(size(A, 2), 1);
newA = zeros(size(A));

for p = 1 : lambda
	%step 2:find maximum A * F 
    max = abs(A(:, 1)' * F);
    index = 1;
    for i = 2 : size(A, 2)
        if(abs(A(:, i)' * F) > max)
            index = i;
            max = abs(A(:, i)' * F);
        end
    end
    %step 3: identity the maximum index
    omiga(p,1) = index;
	%step 4: update omega
    newA(:, omiga(p, 1)) = A(:, omiga(p, 1));
	%step 5: approximate F 
	%step 6: update F
    alpha = newA \ B;
    F = B - newA * alpha;   
end

end