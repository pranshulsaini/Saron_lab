trl=EEG.trials;
code='31_444_2';
for i=1:trl
EEG.epoch(i).eventtype(1)={code};
EEG.event(2*i-1).type=code;
EEG.urevent(i).type=code;
end
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,1);
EEG = eeg_checkset(EEG, 'eventconsistency');