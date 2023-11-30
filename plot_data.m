clc; clear;

load("Data/dp_hjb_1dpoly_n_200_N_300000_epsi0.mat");
%load("Data/dp_hjb_lqr_n_200_N_300000_epsi10_july11.mat");

t_steps = double((0:N))*dt;

%% plots states and control 
%{
fig = figure(1);
plot(t_steps, x_lqr, 'Linewidth',2);
xlabel('time');
ylabel('state - lqr');
saveas(fig,'dp_lqr/x_lqr.png');

fig = figure(2);
plot(t_steps, x_dp, 'Linewidth',2);
xlabel('time');
ylabel('state - DP');
saveas(fig,'dp_lqr/x_dp.png');

fig = figure(3);
plot(t_steps(1:end-1), u_dp, 'Linewidth',2);
xlabel('time');
ylabel('control - DP');
saveas(fig,'dp_lqr/u_dp.png');

fig = figure(4);
plot(t_steps(1:end-1), u_lqr, 'Linewidth',2);
xlabel('time');
ylabel('control - lqr');
saveas(fig,'dp_lqr/u_lqr.png');

fig = figure(5);
plot(t_steps(1:end-1), u_lqr - u_dp, 'Linewidth',2);
xlabel('time');
ylabel('difference in u');
saveas(fig,'dp_lqr/difference_u.png');

fig = figure(6);
plot(t_steps, x_dp -x_lqr, 'Linewidth',2);
xlabel('time');
ylabel('difference in x');
saveas(fig,'dp_lqr/difference_x.png');

%% plot cost-to-go 

[xx, tt]=meshgrid(X,t_steps);
surf(tt,xx,J);
colorbar;
%}

%% comparing cost-to-go of lqr 
n = length(X);
J_lqr = zeros(N+1,n);

for t = 1:N+1
    
    J_lqr(t,:) = 0.5*P(t)*(X.^2) + q(t);
    
end

%% global control of lqr
u_lqr_global = zeros(N,n);

for t = 1:N
    
   u_lqr_global(t,:) = K(t).*X;
    
end

%% plot cost-to-go snapshot. 
VIDEO_OUTPUT = false;

if VIDEO_OUTPUT
    myVideo = VideoWriter('/home/naveed/Documents/Dynamic_programming/cost_to_go_comp_dp_hjb_lqr_n_200_N_300000_epsi10'); %open video file
    myVideo.FrameRate = 10;
    open(myVideo)
end

frame_frequency = 1000;
font_size = 16;
for t=N+1:-frame_frequency:1
%for t=299900:-frame_frequency:1
%for t=245000
    t
    
    fig = figure(1);
    plot(X, J(t,:),'r','LineWidth',3, 'DisplayName', 'HJB-FD');
    hold on;
    %plot(X, J_lqr(t,:),'b--','LineWidth',3, 'DisplayName', 'LQR');
    ax = gca;
    ax.FontSize = 16; 
    xlabel('State','FontSize',font_size);
    ylabel('Cost-to-go','FontSize',font_size);
    %title(['t = ', num2str(t)]);
    h = legend();
    set(h,'FontSize',font_size);
    set(fig,'Units','inches');
    screenposition = get(fig,'Position');
    set(fig,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    %print -dpdf -painters '/home/naveed/Dropbox/Research/Manuscripts/TAC22/plots/ctg_e0.pdf'
    pause(0.001);
    hold off;
    
   
    
    if t ~= N+1
        figure(2);
        plot(X, u_global(t,:),'r','LineWidth',2, 'DisplayName', 'HJB-FD');
        hold on;
        %plot(X, u_lqr_global(t,:),'b--','LineWidth',2, 'DisplayName', 'LQR');
        xlabel('X');
        ylabel('Control');
        title(['t = ', num2str(t)]);
        legend();
        pause(0.01);
        hold off;
    end
    
    
if VIDEO_OUTPUT
%Video output
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);

end
     
end

if VIDEO_OUTPUT
    close(myVideo);
end
