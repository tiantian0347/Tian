clear all

for n = 10 : 10 : 100
    A = (-0.5)*triu(ones(n));
    A = 1.5*eye(n) + A;
    
    fprintf('n=%d, cond(A)=%.4e\n',n,cond(A));
    
end
