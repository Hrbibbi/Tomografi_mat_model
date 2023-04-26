%%
clc, clear, close all
%%
factor = 100;
N = 5000/factor;
load("testImage.mat") % loads array "im"
im_downsample = downsample(im,factor);

figure
imshow(im_downsample,[0 max(max(im_downsample))])
colorbar
%%% noise level always 0.001
% %% varying angles rays = 1*d
% %n_angles = [42,43]; % noise_level = 0
% n_angles = [40,45,50,60,90]; % noise_level = 0.001
% %n_angles = [10,20,30,40,50,100];
% d = sqrt(2)*N;
% p = round(d);
% noise_level = 0.001;
% 
% M = length(n_angles);
% 
% figure
% ha = tight_subplot(1,M,[.01 .03],[.1 .01],[.06 .01]);
% 
% for i = 1:M
%     angles = linspace(0,179,n_angles(i));
%     [im_sol,cond_AtA] = call_tomo(im_downsample,N,angles,p,d,noise_level);
%     disp(n_angles(i)*p)
% 
%     axes(ha(i));
%     if isfinite(max(max(im_downsample)))
%         high = max(max(im_downsample));
%     else
%         high = max(max(im_sol));
%     end
%     disp(high)
%     imshow(im_sol,[0 high])
%     colorbar
%     text(0.5, 0.1, sprintf('cond = %.2e', cond_AtA), ...
%     'HorizontalAlignment', 'center', 'Units', 'normalized', ...
%     'Color', [1 1 0]); % set color to bright yellow
%     if i == 1
%         ylabel(sprintf('p = %d', p), 'Rotation', 0, 'HorizontalAlignment', 'right', 'FontSize', 12, 'FontWeight', 'bold');
%     end
%     title(sprintf('angles = %d', n_angles(i)), 'FontSize', 12)
% end

%% varying angles and rays
% \Delta d is the distance between the rays
% \Delta L is the side-length of pixels
scalers = [1,2,3]; % add 0.5
anglerow = [40,45,50,60,70];
n_angles = round([
    anglerow;
    anglerow/2;
    anglerow/3
    ]);
% n_angles = [
%     40,45,50,60,70;
%     35,37,38,50,60;
%     25,30,35,40,50;
%     ]; % noise_level = 0.001

d = sqrt(2)*N;
p = round(d*scalers);
noise_level = 0.001;

M = length(n_angles);
K = length(p);

figure
ha = tight_subplot(K,M,[.01 .03],[.1 .01],[.07 .01]);

count = 0;
for j = 1:K
    for i = 1:M
        angles = linspace(0,179,n_angles(j,i));
        [im_sol,cond_AtA] = call_tomo(im_downsample,N,angles,p(j),d,noise_level);
        disp(n_angles(i)*p)
        
        count = count + 1;
        axes(ha(count));
        if isfinite(max(max(im_downsample)))
            high = max(max(im_downsample));
        else
            high = max(max(im_sol));
        end
        disp(high)
        imshow(im_sol,[0 high])
        colorbar
        text(0.5, 0.1, sprintf('cond = %.2e', cond_AtA), ...
        'HorizontalAlignment', 'center', 'Units', 'normalized', ...
        'Color', [1 1 0]); % set color to bright yellow
        if i == 1
            ylabel(sprintf('p = %d\n(%cd = %cL/%d)', p(j), 8710, 8710, scalers(j)), 'Rotation', 0, 'HorizontalAlignment', 'right', 'FontSize', 12, 'FontWeight', 'bold');
        end
%         title_str = strcat('nA = ', num2str(n_angles(j,i)), ', prod = ', num2str(p*n_angles(j,i)));
%         title(title_str, 'FontSize', 12)
        title(sprintf('nA = %d,     p*nA = %d', n_angles(j,i), p(j)*n_angles(j,i)), 'FontSize', 12)
%         title(sprintf('nA = %d', n_angles(j,i)), 'FontSize', 12)
    end
end

%% varying angles and rays
% \Delta d is the distance between the rays
% \Delta L is the side-length of pixels
scalers = [0.3,0.6,1];
anglerow = [45,50,60,70,180];
n_angles = [];
for i = 1:length(scalers)
    n_angles = [n_angles; round(anglerow ./ scalers(i))];
end
% n_angles = [
%     40,45,50,60,70;
%     35,37,38,50,60;
%     25,30,35,40,50;
%     ]; % noise_level = 0.001

d = sqrt(2)*N;
p = round(d*scalers);
noise_level = 0.001;

M = length(n_angles);
K = length(p);

figure
ha = tight_subplot(K,M,[.01 .03],[.1 .01],[.09 .01]);
count = 0;
for j = 1:K
    if scalers(j) == 1
        y_str = sprintf('$p = $%d\n$\\Delta d = \\Delta L$', p(j));
    elseif scalers(j) > 1
        y_str = sprintf('$p = $%d\n$\\Delta d = \\Delta L$/%d', p(j), scalers(j));
    else
        y_str = sprintf('$p = $%d\n$\\Delta d =$%d$\\Delta L$', p(j), round(1/scalers(j)));
    end
    for i = 1:M
        angles = linspace(0,179,n_angles(j,i));
        [im_sol,cond_AtA] = call_tomo(im_downsample,N,angles,p(j),d,noise_level);
        disp(n_angles(i)*p)
        
        count = count + 1;
        axes(ha(count));
        if isfinite(max(max(im_downsample)))
            high = max(max(im_downsample));
        else
            high = max(max(im_sol));
        end
        disp(high)
        imshow(im_sol,[0 high])
        colorbar
%         text(0.5, -0.05, sprintf('cond = %.2e', cond_AtA), ...
%         'HorizontalAlignment', 'center', 'Units', 'normalized', ...
%         'Color', [1 1 0]); % set color to bright yellow
        text(0.5, -0.05, sprintf('cond = %.2e', cond_AtA), ...
                'HorizontalAlignment', 'center', 'Units', 'normalized');
%         xlabel(sprintf('cond = %.2e', cond_AtA))
        if i == 1
            ylabel(y_str,'Interpreter', 'latex', 'Rotation', 0, 'HorizontalAlignment', 'right', 'FontSize', 12, 'FontWeight', 'bold');
        end
%         title_str = strcat('nA = ', num2str(n_angles(j,i)), ', prod = ', num2str(p*n_angles(j,i)));
%         title(title_str, 'FontSize', 12)
        title(sprintf('nA = %d,     p*nA = %d', n_angles(j,i), p(j)*n_angles(j,i)), 'FontSize', 12, 'FontWeight', 'normal')
%         title(sprintf('nA = %d', n_angles(j,i)), 'FontSize', 12)
    end
end
%%
exportgraphics(gcf,'figure1.png','Resolution',300,'BackgroundColor','none','ContentType','vector')
%% regularization
% \Delta d is the distance between the rays
% \Delta L is the side-length of pixels
alphas = [0,1,10];
n_angles = [15,30,45,60,180];

d = sqrt(2)*N;
p = round(d);
noise_level = 0.001;

M = length(n_angles);
K = length(alphas);

figure
ha = tight_subplot(K,M,[.01 .03],[.1 .01],[.09 .01]);
count = 0;
for j = 1:K
    for i = 1:M
        angles = linspace(0,179,n_angles(i));
        [im_sol,cond_AtA] = call_tomo_reg(im_downsample,N,angles,p,d,noise_level,alphas(j));
        disp(n_angles(i)*p)
        
        count = count + 1;
        axes(ha(count));
        if isfinite(max(max(im_downsample)))
            high = max(max(im_downsample));
        else
            high = max(max(im_sol));
        end
        disp(high)
        imshow(im_sol,[0 high])
        colorbar
%         text(0.5, -0.05, sprintf('cond = %.2e', cond_AtA), ...
%         'HorizontalAlignment', 'center', 'Units', 'normalized', ...
%         'Color', [1 1 0]); % set color to bright yellow
        text(0.5, -0.05, sprintf('cond = %.2e', cond_AtA), ...
                'HorizontalAlignment', 'center', 'Units', 'normalized');
%         xlabel(sprintf('cond = %.2e', cond_AtA))
        y_str = sprintf('$\\alpha = $%.4f', alphas(j));
        if i == 1
            ylabel(y_str,'Interpreter', 'latex', 'Rotation', 0, 'HorizontalAlignment', 'right', 'FontSize', 12, 'FontWeight', 'bold');
        end
        title(sprintf('nA = %d,     p*nA = %d', n_angles(i), p*n_angles(i)), 'FontSize', 12, 'FontWeight', 'normal')
    end
end

%% Timing using artificial data
clc, clear, close all
%%
%factors = 1000;
N = [20];
% factors = 5000 ./ sizes;
% N = 5000 ./ factors;
angles = linspace(0,179,50);
K = 1:10000;

% load("testImage.mat") % loads array "im"
for i = 1:length(N)
    im_AI = create_data(N(i),1);
    
%     figure
%     imshow(im_downsample,[0 max(max(im_downsample))])
%     colorbar
    
    d = sqrt(2)*N(i);
    p = round(d);
    [T_Gauss(i), T_AIR(i),e_Gauss(i),e_AIR(i),T_BS(i),e_BS(i)] = time_algorithms(im_AI,N(i),angles,p,d,K);
end

% for size = 20 we have T_Gauss = 133.2735 and T_AIR = 0.0058
figure
semilogy(N .^ 2, T_AIR,'b*')
hold on
semilogy(N .^ 2, T_Gauss, 'g*')
semilogy(N .^ 2, T_BS, 'r*')
legend('AIR','Gauss','BS')

%% Timing
clc, clear, close all
%%
%factors = 1000;
sizes = [2,4,5,8,10,20];
factors = 5000 ./ sizes;
N = 5000 ./ factors;
angles = linspace(0,179,50);
K = 1:50;

load("testImage.mat") % loads array "im"
for i = 1:length(N)
    im_downsample = downsample(im,factors(i));
    
%     figure
%     imshow(im_downsample,[0 max(max(im_downsample))])
%     colorbar
    
    d = sqrt(2)*N(i);
    p = round(d);
    [T_Gauss(i), T_AIR(i),e_Gauss(i),e_AIR(i),T_BS(i)] = time_algorithms(im_downsample,N(i),angles,p,d,K);
end

% for size = 20 we have T_Gauss = 133.2735 and T_AIR = 0.0058
figure
semilogy(sizes .^ 2, T_AIR,'bo')
hold on
semilogy(sizes .^ 2, T_Gauss, 'go')
semilogy(sizes .^ 2, T_BS, 'ro')
