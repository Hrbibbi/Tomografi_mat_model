function [T_Gauss, T_AIR,e_Gauss,e_AIR,T_BS,e_BS] = time_algorithms(im,N,angles,p,d,K)
[A,~,~,~,~] = paralleltomo(N,angles,p,d);
x = reshape(im,[numel(im),1]);
b = A*x;
b = b + 0.5*max(b)*randn(size(b));
sys_mat = A'*A;
RHS = (A'*b);
disp('Precomutation complete')

%%% Time Gaussian elimination
% tic
% G = rref(sys_mat);
% x_Gauss = G \ RHS;
% T_Gauss = toc;
% %%%
% disp('Gauss complete')
% e_Gauss = norm(x-x_Gauss,2);
T_Gauss = 0;
e_Gauss = 0;

%%% Time backslash MATLAB
tic
x_BS = A\b;
T_BS = toc;
%%%


e_BS = norm(x-x_BS,2);
disp(e_BS)

%%% Find iteration number for Kaczmarz
[x_AIR,~] = kaczmarz(sys_mat,RHS,K);
for i = 1:length(K)
    if norm(x-x_AIR(:,i),2) < e_BS
        disp('hej')
        disp(i)
        stop_idx = i;
        break
    end
    if i == length(K)
        disp('Fail')
        return
    end
end
disp('K_stop = ' + K(stop_idx))
% [x_AIR,~] = kaczmarz(sys_mat,RHS,K);
% for i = 1:length(K)
%     if norm(x-x_AIR(:,i),2) < e_Gauss
%         disp('hej')
%         disp(i)
%         stop_idx = i;
%         break
%     end
%     if i == length(K)
%         disp('Fail')
%         return
%     end
% end
%disp('K_stop = ' + K(stop_idx))

%%% Time Kaczmarz
%disp(K(stop_idx))
tic
[x_AIR,~] = kaczmarz(sys_mat,RHS,K(stop_idx));
T_AIR = toc;
%%%

e_AIR = norm(x-x_AIR,2);