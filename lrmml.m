function[Wtplus, Rtplus]=lrmml(Xtrain, Ytrain, lambda, lambda1, lambda2)
%Code based on LSML and DM2L multi-label learning algorithms.
%Requires poblano_toolbox_1.1
[n, m]=size(Ytrain);
[~, d]=size(Xtrain);
X=cell(m, 1);
for i=1:m
    tempindex=find(Ytrain(:, i)>0);
    X{i}=Xtrain(tempindex, :);    
end
J = Ytrn ~= 0;
XTX = Xtrain' * Xtrain;
XTY = Xtrain' * Ytrain;
num_dim = d;
rho = 2^1;
Wt = (XTX + rho*eye(num_dim)) \ (XTY);
Rt = eye(m, m); 
Flag=true;
iteration=1;
while Flag && iteration <25
    [Wtplus, Rtplus]=optsurrogate(Xtrain, X, Ytrain, J, Wt, Rt, lambda, lambda1, lambda2);
    if norm(Wtplus-Wt, 'fro')/norm(Wt, 'fro') <10^-3
        Flag=false;
    else
       iteration=iteration+1; 
       Wt=Wtplus;
       Rt=Rtplus;
    end  
end
end

