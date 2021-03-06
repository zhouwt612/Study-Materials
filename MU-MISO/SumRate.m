function C = SumRate(H,P,sigma2)
[K,M] = size(H);
c = zeros(1,K);
for idx1 = 1:1:K
    ds = abs(H(idx1,:)*P(:,idx1))^2;
    int = 0;
    for idx2 = 1:1:K
        if idx2 ~= idx1
            int = int + abs(H(idx1,:)*P(:,idx2))^2;
        end
    end
    sinr_k = ds/(sigma2+int);
    c(idx1) = log2(1+sinr_k);
end
C = sum(c);