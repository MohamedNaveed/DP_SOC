% plot cost of DP and LQR

epsilon = [0, 0.2, 0.4, 0.6, 0.8, 1.0];

DP_cost = [5.9486138465, 6.1935342295, 7.3588599796, 9.1539483023, ...
            11.6905510559, 14.3681037985];
      
DP_var = [3.15544362088405e-30, 0.2424658113, 1.2561649144, 4.4554127424,...
            11.4203340717, 20.1542351755];

DP_std = sqrt(DP_var);

LQR_cost = [5.8172454104, 6.1577496237, 7.5294154444, 9.3288781859,...
            11.4086018647, 15.6054890876];
LQR_var = [7.88860905221012e-31, 0.264352348, 2.4701954279, 4.6960852161,...
            8.3623044541, 27.4549519426];
LQR_std = sqrt(LQR_var);

fontsize = 14;

fig = figure;
curve1 = DP_cost + DP_std;
curve2 = DP_cost - DP_std;
x2 = [epsilon, fliplr(epsilon)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r', 'FaceAlpha',0.1,'LineStyle','none');
hold on;
plot(epsilon, curve1,'r');
plot(epsilon, curve2,'r');

curve1 = LQR_cost + LQR_std;
curve2 = LQR_cost - LQR_std;
x2 = [epsilon, fliplr(epsilon)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b','FaceAlpha',0.1,'LineStyle','none');
hold on;
plot(epsilon, curve1,'b');
plot(epsilon, curve2,'b');
plot(epsilon, DP_cost,'r', 'LineWidth',2, 'DisplayName', 'DP');
plot(epsilon, LQR_cost,'--b', 'LineWidth',2, 'DisplayName', 'LQR');

xlabel('epsilon', 'Fontsize', fontsize);
ylabel('Cost','Fontsize', fontsize);
f=get(gca,'Children')
h = legend([f(1),f(2)],'LQR','DP');
set(h,'FontSize',fontsize);

set(fig,'Units','inches');
fig.Position = [80,80,4.5,4.5];
screenposition = get(fig,'Position');
set(fig,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters '/home/naveed/Documents/Dynamic_programming/DP_LQR_cost.pdf'





    





