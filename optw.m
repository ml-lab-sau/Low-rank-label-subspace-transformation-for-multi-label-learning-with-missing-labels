function [f,g] = optw(x,Data)
[d, m] = size(Data.W);
Wplus=reshape(x, d, m);
E = ones(size(Data.Ytrain));
HingeL = max((E - (Data.Ytrain* Data.R) .* (Data.Xtrain * Wplus)) , 0);
term1 = 1/2 * norm(Data.J.*(HingeL.^2),1);
term1R = Data.lambda1 * 1/2 * norm(Data.J .* (Data.Ytrain * Data.R - Data.Ytrain), 'fro')^2;
Sabs=0;
for i=1:m
  S=svd(Data.X{i}*Wplus, 'econ');
  Sabs=Sabs+sum(abs(S));
end
term2=Data.lambda*Sabs;
Normgradient=subgradient_nuclearnorm(Data.Xtrain*Data.W);
term3=Data.lambda*trace(Wplus'*Data.Xtrain'*Normgradient);
L = diag(sum(Data.R,2)) - Data.R;
term4 = Data.lambda2 * trace(Data.W*L*Data.W');
f = term1+term2-term3+term4;
gterm1=Data.Xtrain' * (HingeL .* (-(Data.Ytrain * Data.R) .* Data.J));
gterm2=zeros(d, m);
for i=1:m
   gterm2=gterm2+Data.X{i}'*subgradient_nuclearnorm(Data.X{i}*Wplus);  
end
gterm2=Data.lambda*gterm2;
gterm3=Data.lambda*Data.Xtrain'*subgradient_nuclearnorm(Data.Xtrain*Data.W);
gterm4=Data.lambda2 * Data.W *(L + L');
g=gterm1+gterm2-gterm3+gterm4;
g=reshape(g, d*m, 1);
end