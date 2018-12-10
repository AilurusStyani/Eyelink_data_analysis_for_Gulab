function bar_error_plot(figurei,y1,y2,style,name)

y1=y1(:);
y2=y2(:);
if ishandle(figurei); clf;close (figurei); end
figure(figurei);
set(figure(figurei),'Name',name,'color','white'); 

b_bar(1)=bar(1,nanmean(y1));
hold on;

set(b_bar(1),'FaceColor','b');
b_bar(2)=bar(2,nanmean(y2));
set(b_bar(2),'FaceColor','r');
hold on;

if style(2)==1
    std_y1=nanstd(y1);
    std_y2=nanstd(y2);
    errorbar(1,nanmean(y1),std_y1,'b','LineWidth',4,'CapSize',10);
    errorbar(2,nanmean(y2),std_y2,'r','LineWidth',4,'CapSize',10);

end

if style(1)==1
%     if length(y1)>20000
%         plot(ones(50000),y1(1:50000),'ko')
%         plot(ones(length(y1)-50000),y1(50000:end),'ko');
%     else
        plot(ones(length(y1),1),y1,'ko');
        hold on;
%     end
%     if length(y2)>50000
%         plot(ones(50000),y2(1:50000),'ko')
%         plot(ones(length(y2)-50000),y2(50000:end),'ko');
%     else
        plot(ones(length(y2),1)+1,y2,'ko');
        hold on;
%     end
end


ylim([0 inf]);
% ylim([0 2500]);
title(name);
end