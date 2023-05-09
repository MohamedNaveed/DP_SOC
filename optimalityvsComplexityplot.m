fig = figure(1);
hold on;
xlabel('Complexity','FontSize', 12);
ylabel('Practical Optimality','FontSize', 12);
xlim([0,1]);
ylim([0,1]);
set(gca,'XTick',[],'YTick',[]);
x = [0.63,0.40,0.03];
y = [0.9,0.75,0.58];
txt = {'MPC',['Linear Feedback',char(10), 'with replanning'],'Linear Feedback'};
%plot(0.99,0.82,'ro','Markersize',10,'MarkerFaceColor','r');
plot(0.75,0.9,'bo','Markersize',10,'MarkerFaceColor','b');
plot(0.45,0.85,'go','Markersize',10,'MarkerFaceColor','g');
plot(0.05,0.65,'ko','Markersize',10,'MarkerFaceColor','r');

text(x,y,txt,'FontSize',14);

grid on 
pbaspect([1.5 1 1])
set(fig,'Units','inches');
screenposition = get(fig,'Position');
set(fig,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
grid on 
%print -dpdf -painters '/home/naveed/Documents/WAFR_presentation/OptimalityvsTractability.pdf' ;

%hold off;