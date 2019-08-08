% header_list = ["x_data", "y_data", "error_0", "error_1"]

openfig('IrO2round1.fig');

x_data = get(get(gca, 'Children'), 'XData');
y_data = get(get(gca, 'Children'), 'YData');

x_data = x_data(2);
y_data = y_data(2);

x_data = x_data{1,1};
y_data = y_data{1,1};


errorbar = get(get(gcf, 'Children'), 'Children');
errorbar = errorbar(2);

error_0 = get(errorbar, 'YNegativeDelta');
error_1 = get(errorbar, 'YPositiveDelta');


tmp0 = [x_data; y_data; error_0; error_1]';

csvwrite('IrO2round1.csv', tmp0)
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

openfig('IrO2round2.fig');

x_data = get(get(gca, 'Children'), 'XData');
y_data = get(get(gca, 'Children'), 'YData');

x_data = x_data(2);
y_data = y_data(2);

x_data = x_data{1,1};
y_data = y_data{1,1};


errorbar = get(get(gcf, 'Children'), 'Children');
errorbar = errorbar(2);

error_0 = get(errorbar, 'YNegativeDelta');
error_1 = get(errorbar, 'YPositiveDelta');


tmp0 = [x_data; y_data; error_0; error_1]';

csvwrite('IrO2round2.csv', tmp0)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

openfig('IrO2round3.fig');

x_data = get(get(gca, 'Children'), 'XData');
y_data = get(get(gca, 'Children'), 'YData');


tmp0 = get(get(gcf, 'Children'), 'Children');
error_0 = get(tmp0, 'YNegativeDelta');
error_1 = get(tmp0, 'YPositiveDelta');


tmp0 = [x_data; y_data; error_0; error_1]';

csvwrite('IrO2round3.csv', tmp0)
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

