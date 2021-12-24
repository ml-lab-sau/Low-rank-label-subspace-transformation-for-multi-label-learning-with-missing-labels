function [f,g] = optr(x,Data)
[~, m] = size(Data.R);
Rplus=reshape(x, m, m);
E = ones(size(Data.Ytrain));
HingeL = max((E - (Data.Ytrain * Rplus).* (Data.Xtrain * Data.W)) , 0);
term1 = 1/2 * norm(Data.J.*(HingeL.^2),1);
term1R = Data.lambda1 * 1/2 * norm(Data.J.*(Data.Ytrain * Rplus - Data.Ytrain),'fro')^2;
Sabs=0;
for i=1:m
  S=svd(Data.X{i}*Data.W, 'econ');
  Sabs=Sabs+sum(abs(S));
end
term2=Data.lambda*Sabs;
Normgradient=subgradient_nuclearnorm(Data.Xtrain*Data.W);
term3=Data.lambda*trace(Data.W'*Data.Xtrain'*Normgradient);
L = diag(sum(Data.R,2)) - Data.R;
term4 = Data.lambda2 * trace(Data.W*L*Data.W');
f = term1+term1R+term2-term3+term4;
gterm1=Data.Ytrain' * (HingeL .* (-(Data.Xtrain * Data.W) .* Data.J));
gterm1R = Data.lambda1 * Data.Ytrain' * ((Data.Ytrain * Rplus - Data.Ytrain) .* Data.J);
g=gterm1+gterm1R;
g=reshape(g, m*m, 1);
end