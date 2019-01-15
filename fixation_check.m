function dataout = fixation_check(datain,degree)
global SCREEN
distance = sqrt((datain(:,2) - SCREEN.screenHeight / 2).^2 + (datain(:,3) - SCREEN.screenWidth / 2).^2);
if SCREEN.screenWidth / SCREEN.screenWidth_real ~= SCREEN.screenHeight / SCREEN.screenHeight_real
    error('��ǰ��������������ʾ�������������Ĳ����������ʾ���Ƿ�Ϊ����')
end
valid_distance = SCREEN.view_distance * tand(degree) / SCREEN.screenWidth_real * SCREEN.screenWidth;
datain(distance-valid_distance>0,2:end)=NaN;
dataout = datain;