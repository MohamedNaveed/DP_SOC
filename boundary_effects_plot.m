% plot boundary effect. 
clc;clear;close all;
%load("Data/dp_hjb_1dcos_n_200_N_300000_epsi0_domain_pi.mat");
load("Data/dp_hjb_lqr_n_200_N_300000_epsi10_july11.mat");

%%
t_steps = double((0:N))*dt;

x0 = 1;
n_samples = 10;
x = zeros(N+1,n_samples);

epsilon = 10;

for i = 1:n_samples
   
    x(1,i) = x0;
    
    for t = 1:N
        x_idx = find_nearest(X,x(t,i));
        u = u_global(t,x_idx);

        x(t+1,i) = model(x(t,i), u, dt, epsilon);

    end
end

%%

fig = figure(1);
hold on;
t_idx = 1:100:N+1;
plot(t_steps(t_idx), x(t_idx,1:10), 'LineWidth',2, 'HandleVisibility','off');
%plot(t_steps, x(:,2), '--', 'LineWidth',2, 'DisplayName', 'epsilon = 0.2');
%plot(t_steps, x(:,3), '--', 'LineWidth',2, 'DisplayName', 'epsilon = 0.8');
xlabel('time');
ylabel('state');
plot(t_steps(t_idx), 2*ones(length(t_idx),1),'k','LineWidth', 2, 'DisplayName', 'Domain boundary');
plot(t_steps(t_idx), -2*ones(length(t_idx),1),'k','LineWidth', 2, 'HandleVisibility','off');
ylim([-3,3]);
grid on;
%h = legend();
font_size = 16;
ax = gca;
ax.FontSize = font_size; 
%set(h,'FontSize',font_size);
set(fig,'Units','inches');
screenposition = get(fig,'Position');
set(fig,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
%print -dpdf -painters '/home/naveed/Dropbox/Research/Manuscripts/TAC22/plots/traj_epsi10.pdf'



















function [x_new] = model(x,u,dt, epsilon)

    x_dot = x + u ;
    x_new = x + x_dot*dt + epsilon*normrnd(0,1)*sqrt(dt);
end

function [idx] = find_nearest(inp_arr, value)

[min_val, idx] = min(abs(inp_arr - value));

end