neurons = 5;
window = 1:10000;

figure
plot(neural_data_temp_M1(neurons,window)'); hold on
speed = sqrt(x_vel(window).^2+y_vel(window).^2);
normz = max(max(neural_data_temp_M1(neurons,window)))/max(speed);
plot(speed*normz,'LineWidth',3); hold off


figure
hist(kin.pos(:,1))
figure
hist(ts)