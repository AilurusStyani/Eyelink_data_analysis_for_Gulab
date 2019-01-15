%far,close,distance,choice,time
clear all
clc
clf
load('D:\data0105\heading\DQS\microsaccade_direction2.csv');

delect = find(microsaccade_direction2(:,1)==0&microsaccade_direction2(:,2)==0);
microsaccade_direction2(delect,:)=[];
microsaccade_direction2(:,6)= (microsaccade_direction2(:,1)+microsaccade_direction2(:,2));
far_porpotion = microsaccade_direction2(:,1)./microsaccade_direction2(:,6);

index = find(microsaccade_direction2(:,4)==2);
far_choice = far_porpotion(index);
index = find(microsaccade_direction2(:,4)==1);
close_choice = far_porpotion(index);


microsaccade_dire = microsaccade_direction2(:,1)==1;
[tbl,chi2,p,labels] = crosstab(microsaccade_dire, microsaccade_direction2(:,4)) 
set(0,'defaultfigurecolor','w');
% for j = 1:length(far_choice)
% figure(1);
% hold on
% plot(far_choice(j),0.0005*j,'.r');
% % set(gca,'YLim',[0 2]);
% end
% for k = 1:length(close_choice)
% figure(1);
% hold on
% plot(close_choice(k),0.0005*k,'.b');
% % set(gca,'YLim',[0 2]);
% end
% ttest2(far_choice,close_choice)
% set(0,'defaultfigurecolor','w');
% figure(2);
% hold on
% far_choice_porpotion = tabulate(far_choice);
% bar(far_choice_porpotion(:,1),far_choice_porpotion(:,3),'b');
% figure(3);
% close_choice_porpotion = tabulate(close_choice);
% bar((close_choice_porpotion(:,1)),close_choice_porpotion(:,3),'r');
% 
% mean(far_choice)
% mean(close_choice)

%    plot(microSaccade_index(:,1), ones(length(microSaccade_index(:,1)),1)*5+j*0.05,'.k');
