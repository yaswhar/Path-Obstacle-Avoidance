%% Yashar Zafari - 99106209
%% Question 1 & 2
Xs=[1;10];
Xf=[22;12];
eta=1;
alpha=1;
eps=0.1;
p_o=2;
B=[[4 12 1]' [10 18 1]' [12 16 1]' [6 10 1]'...
   [9 10 2]' [12 4 2]' [6 4 2]'...
   [11 12 3]' [13 16 3]' [17 16 3]' [19 12 3]' [15 8 3]'];
obs=polyshape;
for i=1:3
    obs(i)=polyshape(B(1:2,B(3,:)==i)');
end
P=Path_generator(Xs,Xf,eta,B,alpha,eps,p_o);
plot(obs)
hold on
plot(P(1,:),P(2,:))
plot(P(1,[1 end]),P(2,[1 end]),'Marker','x','LineStyle','none','Color','k')
axis equal
grid on
title('Planned Path')
labels={'Start Point' 'Final Point'};
text(P(1,[1 end]),P(2,[1 end]),labels,"VerticalAlignment","cap")
%% Question 3
B=[[4 12 1]' [10 18 1]' [12 16 1]' [6 10 1]'...
   [9 10 2]' [12 4 2]' [6 4 2]'...
   [11 12 3]' [13 16 3]' [17 16 3]' [19 12 3]' [15 8 3]'...
   [17.5 10 4]' [17.5 8 4]'];
P=Path_generator(Xs,Xf,eta,B,alpha,eps,p_o);
figure
plot(obs)
hold on
plot([17.5 17.5],[10 8]);
plot(P(1,:),P(2,:))
title('Stuck in Local Minimum')
plot(P(1,[1 end]),P(2,[1 end]),'Marker','x','LineStyle','none','Color','k')
axis equal
grid on
labels={'Start Point' 'Local Minimum'};
text(P(1,[1 end]),P(2,[1 end]),labels,"VerticalAlignment","cap")
%% Question 4
B=[[4 12 1]' [10 18 1]' [12 16 1]' [6 10 1]'...
   [9 10 2]' [12 4 2]' [6 4 2]'...
   [11 12 3]' [13 16 3]' [17 16 3]' [19 12 3]' [15 8 3]'...
   [17.5 10 4]' [17.5 8 4]'];
P=Path_generator_esc(Xs,Xf,eta,B,alpha,eps,p_o);
figure
plot(obs)
hold on
plot([17.5 17.5],[10 8]);
plot(P(1,:),P(2,:))
title('Escaping the Local Minimum')
plot(P(1,[1 end]),P(2,[1 end]),'Marker','x','LineStyle','none','Color','k')
axis equal
grid on
labels={'Start Point' 'ّFinal Point'};
text(P(1,[1 end]),P(2,[1 end]),labels,"VerticalAlignment","cap")