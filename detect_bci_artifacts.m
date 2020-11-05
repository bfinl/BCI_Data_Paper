function [bad_chan, artifact_bool] = detect_bci_artifacts(BCI)
%{
function detects artifacts in raw EEG data
Bad channels are defined as those with variance > 5z compared to other
channels
artifacts are defined as trials with any values over 100 microvoltz
requirements: Fieldtrip-201904019
input: BCI-BCI data structure
output: bad_chan - vector of noisy channels
        artifact_bool - vector 1 if trial contains artifact 0 if trial does
        not contain artifact
%}

%get percentage of non-artifact trials
%add obj information
bad_chan = [];
artifact_bool = [];
BCI_runs = unique([BCI.TrialData(:).runnumber]);
for run = 1:length(BCI_runs)
    run_ind = find([BCI.TrialData(:).runnumber] == BCI_runs(run));
    trials = BCI.data(run_ind);
    artifact_found = zeros(length(trials),1);
    
    %look for bad electrodes
    artifact_elec = 1;
    bad_chan_run = [];
    while artifact_elec
        trials_full = cat(2,trials{:});
        trials_full(bad_chan_run,:) = 0;
        trials_full_filt = ft_preproc_bandpassfilter(trials_full,BCI.SRATE,[8 30]);
        elec_var = var(trials_full_filt,[],2);
        elec_z = zscore(elec_var);
        tmp_artifact = find(elec_z > 5);
        if isempty(tmp_artifact)
            artifact_elec = 0;
        else
            bad_chan_run = [bad_chan_run;tmp_artifact];
        end
    end
    
    for trial = 1:length(trials)
        trial_filt = ft_preproc_bandpassfilter(trials{trial},BCI.SRATE,[8 30]);
        %constrict to feedback exclusively
        trial_filt = trial_filt(:,(4*BCI.SRATE + 1):BCI.TrialData(run_ind(trial)).resultind);
        trial_filt(bad_chan_run,:) = 0;
        artifact_found(trial) = logical(sum(abs(trial_filt(:)) > 100));
    end
    artifact_bool = [artifact_bool;artifact_found];
    bad_chan = [bad_chan;bad_chan_run];
end%for BCI run

bad_chan = unique(bad_chan);
           
end%function