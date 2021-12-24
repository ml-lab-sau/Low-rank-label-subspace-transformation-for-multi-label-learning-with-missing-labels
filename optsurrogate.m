function[Wplus, Rplus]=optsurrogate(Xtrain, X, Ytrain, J, W, R, lambda, lambda1, lambda2)
Data.Xtrain=Xtrain;
Data.X=X;
Data.Ytrain=Ytrain;
Data.J=J;
Data.W=W;
Data.R=R;
Data.lambda=lambda;
Data.lambda1=lambda1;
Data.lambda2=lambda2;
[d, m]=size(W);
[m, m]=size(R);
x0=reshape(W, d*m ,1);
 out=ncg(@(x) optw(x, Data), x0, 'MaxIters', 1, 'RelFuncTol', 1e-5, 'StopTol', 1e-6, ...)
'MaxFuncEvals', 100, 'Display', 'final');
Wplus=reshape(out.X, d, m);
Data.W = Wplus;
x1=reshape(R, m*m ,1);
 out=ncg(@(x) optr(x, Data), x1, 'MaxIters', 1, 'RelFuncTol', 1e-5, 'StopTol', 1e-6, ...)
'MaxFuncEvals', 100, 'Display', 'final');
Rplus=reshape(out.X, m, m);
Data.R = Rplus;
end

