function [Bwmmse,SR] = WMMSEpreandSR(H,Pinit,Etx,u,P,Q,K,itenum)
B = Pinit'*sqrt(Etx/trace(Pinit'*Pinit));
Rn = nointercov(H,B,P,Q,K);
% WSRBF-WMMSE
Ammse = zeros(Q*K,Q*K);
W = zeros(Q*K,Q*K);
Datarate = zeros(1,itenum);
Datarate(1) = SUMrate2(H,B,Rn,u,Q,K);
for index = 1:1:itenum

    % Update Ammse
    for i = 1:1:K
        Ammse(((i-1)*Q+1):i*Q,((i-1)*Q+1):i*Q) = RxMMSE(H(((i-1)*Q+1):i*Q,:),...
            B(:,((i-1)*Q+1):i*Q),Rn(:,:,i));
    end
    % Update W
    for j = 1:1:K
        W(((j-1)*Q+1):j*Q,((j-1)*Q+1):j*Q) = MSEweight(H(((j-1)*Q+1):j*Q,:)...
            ,B(:,((j-1)*Q+1):j*Q),Rn(:,:,j));
    end
    % Update Bmmse
    B = TxWMMSE(H,Ammse,W,Etx);

    % Compute the sum rate
    Rn = nointercov(H,B,P,Q,K);
    Datarate(index+1) = SUMrate2(H,B,Rn,u,Q,K);
    Bwmmse = zeros(P,Q,K);
end
SR = Datarate;
Bwmmse = zeros(P,Q,K);
for idx2 = 1:1:K
    Bk = B(:,((idx2-1)*Q+1):idx2*Q);
    Bwmmse(:,:,idx2) = Bk;
end
% Bwmmse = B;


function A_mmse = RxMMSE(H,B,Rn)
% H (QxP): Channel realization
% B (PxQ): Transmitter filter (Beamformer)
% Rn (QxQ): Noise covariance matrix
A_mmse = B'*H'/(H*B*B'*H' + Rn);

function Bwmmse = TxWMMSE(H,A,W,Etx)
% H (QKxP): Channel realization
% A (QKxQK): MMSE receive filter
% W (QKxQK): MSE weights
% Etx: Power constraint
% B (BxQK): Receive filter
B_tilde = inv(H'*A'*W*A*H + (trace(W*A*A')/Etx)*eye(size(H'*H)))*H'*A'*W;
b = sqrt(Etx/trace(B_tilde*B_tilde'));
Bwmmse = b*B_tilde;

function W = MSEweight(H,B,Rn,u)
% E: MMSE matrix
% H (QxP): Channel realization
% B (PxQ): Transmitter filter (Beamformer)
% Rn (QxQ): Noise covariance matrix
% u: parameter
if nargin < 4
    u = 1;
end
E = inv(eye(size(H*H')) + B'*H'/Rn*H*B);
W = u*inv(E);

function Rvv = nointercov(H,B,P,Q,K)
% Rn: Effective noise covariance matrix at user k
Rvv = zeros(Q,Q,K);
B_ite = zeros(P,Q,K);
H_ite = zeros(Q,P,K);

for j = 1:1:K
    B_ite(:,:,j) = B(:,(j-1)*Q+1:j*Q);
    H_ite(:,:,j) = H((j-1)*Q+1:j*Q,:);
end

for k = 1:1:K
    for i = 1:1:K
        if i ~= k
            Rvv(:,:,k) = Rvv(:,:,k) + H_ite(:,:,k)*B_ite(:,:,i)*B_ite(:,:,i)'*H_ite(:,:,k)';
        else
            Rvv(:,:,k) = Rvv(:,:,k) + zeros(Q);
        end
    end 
end
% A normal case, Rvv = Rvv + sigma2*eye(Q);
% But the noise power is set to 1 in our model.
Rvv = Rvv + eye(Q);

function sumrate = SUMrate2(H,B,Rn,u,Q,K)
% E: MMSE matrix
% H (QxP): Channel realization
% B (PxQ): Transmitter filter (Beamformer)
% Rn (QxQ): Noise covariance matrix
% u: parameter
% k: Number of users
R = zeros(1,K);
for idx = 1:1:K
    Bk = B(:,((idx-1)*Q+1):idx*Q);
    Hk = H(((idx-1)*Q+1:idx*Q),:);
    Ek = eye(size(Hk*Hk')) + Bk'*Hk'/Rn(:,:,idx)*Hk*Bk;
    Rk = log2(det(Ek));
    R(idx) = u(idx)*Rk;
end
sumrate = real(sum(R));