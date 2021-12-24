function[deltaA]=subgradient_nuclearnorm(A)
% An m * n matrix A
[m, n]=size(A);
[U, B, V]=svd(A);
s=length(find(diag(B)>=0.005));
if m <s || n <s
    deltaA=zeros(m, n);
else
    U1=U(:, 1:s);
    V1=V(1:s, :);
    deltaA=U1*V1;
end
end