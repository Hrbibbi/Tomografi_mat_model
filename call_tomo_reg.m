function [im_sol,cond_AtA] = call_tomo(im,N,angles,p,d,noise_level,alpha)
[A,~,~,~,~] = paralleltomo(N,angles,p,d);%,0.01)%,round(sqrt(2)*N),sqrt(2)*N);%,0.1)
x_true = reshape(im,[numel(im),1]);
b = A*x_true;

b_noise = b + noise_level*max(b)*randn(size(b));
x_sol = ( A'*A + alpha*eye(size(A,2)) )\(A'*b_noise);
%x_sol = A\b_noise;
im_sol = reshape(x_sol,[N,N]);

cond_AtA = cond(A'*A + alpha*eye(size(A,2)));