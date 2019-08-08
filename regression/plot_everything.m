

X = [EEG.chanlocs.Y];
Y = [EEG.chanlocs.X];
Z = [EEG.chanlocs.Z];
elec = {EEG.chanlocs.labels};
id_elec = find(strcmp('Cz', elec));

figure
hold on
plot3([EEG.chanlocs.Y], [EEG.chanlocs.X], [EEG.chanlocs.Z], '.')
view([0, 0, 1])
textscatter3(X(id_elec), Y(id_elec), Z(id_elec), elec(id_elec))

% % % for i = 1:63
% % % textscatter3(X(i), Y(i), Z(i), elec(i))
% % % end


