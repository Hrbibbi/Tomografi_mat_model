%%
% X = reshape(1:16,4,4)
% R = reshape(X, 2, 2, 2, 2);
% S = sum(sum(R, 1), 3) * 0.25;
% Y = reshape(S, 2, 2);

% load("testImage.mat")
% % X = reshape(1:64,8,8)
% % Y = downsample(X,2)
% im_downsample = downsample(im,100)
% imshow(im_downsample,[0 max(max(im_downsample))])

% lead: 0.0038
% wood: 0.0034
% steel: 0.0100

%%
% https://se.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w
%%
clc, clear, close all
%%
factor = 100;
N = 5000/factor;
load("testImage.mat") % loads array "im"
im_downsample = downsample(im,factor);
% x_true = reshape(im_downsample,[numel(im_downsample),1]);


figure
imshow(im_downsample,[0 max(max(im_downsample))])
colorbar
%% varying angles
%n_angles = [42,43]; % noise_level = 0
n_angles = [40,45,50,60,180]; % noise_level = 0.001
%n_angles = [10,20,30,40,50,100];
d = sqrt(2)*N;
p = round(d);
noise_level = 0.001;

M = length(n_angles);

figure
% ha = tight_subplot(1,M,[.01 .03],[.1 .01],[.01 .01]);
ha = tight_subplot(1,M,[.01 .03],[.1 .01],[.06 .01]);

for i = 1:M
    angles = linspace(0,179,n_angles(i));
    [im_sol,cond_AtA] = call_tomo(im_downsample,N,angles,p,d,noise_level);
    disp(n_angles(i)*p)

    %subplot(1,M,i)
    %subaxis(1,M,i, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.05);
    axes(ha(i));
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
%     text(0.5, 0.5, sprintf('Cond: %.2e', cond_AtA), ...
%         'HorizontalAlignment', 'center', 'Units', 'normalized');
    if i == 1
        ylabel(sprintf('p = %d', p), 'Rotation', 0, 'HorizontalAlignment', 'right', 'FontSize', 12, 'FontWeight', 'bold');
    end
    title(sprintf('angles = %d', n_angles(i)), 'FontSize', 12)%, 'FontWeight', 'normal');
end
% set(ha(1:N),'XTickLabel',''); set(ha,'YTickLabel','')

%%
% example values of the two parameters
a = [1, 2, 3];
b = [10, 20];

% calculate the number of rows and columns of the subplot array
m = length(b);
n = length(a);

% create a figure
figure;

% loop over all subplots
for i = 1:m*n
    % calculate the row and column indices of the current subplot
    row = ceil(i/n);
    col = mod(i-1,n) + 1;
    
    % create the subplot and plot some data
    subplot(m,n,i);
    plot(rand(1,10)*row + col);
    
    % remove the x and y tick labels
    set(gca, 'XTickLabel', '');
    set(gca, 'YTickLabel', '');
    ylabel('', 'Rotation', 0, 'HorizontalAlignment', 'right')
    % add x and y axis labels to the bottom and leftmost subplots
    if row == m
        xlabel('Parameter 1');
    end
    if col == 1
        ylabel('Parameter 2');
    end
end
%%
% define parameter ranges
param1_range = [1 2 3];
param2_range = [4 5 6];

% create subplot array
figure;
ha = tight_subplot(numel(param2_range), numel(param1_range), 0.04, [0.1 0.05], [0.01 0.01]);

% loop over parameter values and plot
for i = 1:numel(param1_range)
    for j = 1:numel(param2_range)
        % compute data for current parameter values
        x = 0:0.01:1;
        y = sin(2*pi*x*param1_range(i)) .* cos(2*pi*x*param2_range(j));
        
        % select current subplot and plot data
        axes(ha((j-1)*numel(param1_range)+i));
        plot(x,y);
        
        % set axis labels
        if j == 1
            title(sprintf('Parameter 1 = %d', param1_range(i)), 'FontWeight', 'normal');
        end
        if i == 1
            ylabel(sprintf('Parameter 2 = %d', param2_range(j)));
        end
        
        % adjust axis properties
        set(gca, 'YTick', []);
        set(gca, 'XTickLabel', '');
    end
end

% set y-label at the top
set(ha, 'YDir','reverse');
set(ha, 'YTick', fliplr(param2_range));

%%
% define parameter ranges
param1_range = [1 2 3];
param2_range = [4 5 6];

% create subplot array
figure;
ha = tight_subplot(numel(param2_range), numel(param1_range), 0.04, [0.1 0.05], [0.05 0.01]);

% loop over parameter values and plot
for i = 1:numel(param1_range)
    for j = 1:numel(param2_range)
        % compute data for current parameter values
        x = 0:0.01:1;
        y = sin(2*pi*x*param1_range(i)) .* cos(2*pi*x*param2_range(j));
        
        % select current subplot and plot data
        axes(ha((j-1)*numel(param1_range)+i));
        plot(x,y);
        
        % set axis labels
        if j == 1
            title(sprintf('Parameter 1 = %d', param1_range(i)), 'FontWeight', 'normal');
        end
        if i == numel(param1_range)
            ylabel(sprintf('Parameter 2 = %d', param2_range(j)), 'FontWeight', 'normal', 'HorizontalAlignment', 'right');
        end
        
        % adjust axis properties
        set(gca, 'YTick', []);
        set(gca, 'XTickLabel', '');
    end
end

% set y-label at the top
set(ha, 'YDir','reverse');
set(ha, 'YTick', fliplr(param2_range));

%%

% figure
% fig=gcf;
% fig.Units='normalized';
% fig.OuterPosition=[0 0 1 1];

% b_noise = b + 0.001*max(b)*randn(size(b));
% x_sol = A\b_noise;
% im_sol = reshape(x_sol,[N,N]);
% imshow(im_sol,[0 max(max(im_downsample))])
% colorbar

% b_noise = imnoise(b,'gaussian')
% sol = A\b_noise
% im_sol = reshape(sol,[30,30])

% imshow(im_x,[min(min(im_x)) max(max(im_x))])
% figure
% imshow(im_sol,[min(min(im_sol)) max(max(im_sol))])