
[dCOPy, dFy] = meshgrid(-50:50, -10:0.2:10);
hwidth = 29.6 * 2;
dCOPz = hwidth; Mcom = 395;
Fg = (dCOPz * dFy + Mcom) ./ dCOPy;
figure
surf(dFy, dCOPy, Fg)

hold on
Mcom = -395;
Fg = (dCOPz * dFy + Mcom) ./ dCOPy;
surf(dFy, dCOPy, Fg)

%%
xlim([-10, 15])
ylim([-30, 20])
zlim([15, 35])
grid on
set(gca, 'View', [70, 35])
% % % axis equal
xlabel('{\Delta}COPy_{TH-VF} (mm)');
ylabel('{\Delta}Fy_{TH-VF} (N)');
zlabel('Fn_{TH} (N)');