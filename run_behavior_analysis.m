function run_behavior_analysis(data_path)
%{
function collects behavioral data from the BCI files
input:  data_path: were datafiles are stored

output: figures 1 and 2 from data paper
%}
%script setup

path.data = data_path;
subfiles = dir(path.data);
subfiles = {subfiles(~[subfiles.isdir]).name};

%metadata hold information regarding dataset
metadata.subjects = 62;
metadata.sessions = length(subfiles);
md_fields = {'runs','trials','total_time','feedback_time'};%time in seconds

for mdf = 1:length(md_fields)
    metadata.(md_fields{mdf}) = 0;
end

%data matricies to hold collected information
%bci_pvc: subs, sessions, pdm
bci_pvc = NaN(metadata.subjects,11,3);

%percent of artifact trials
artifact_percentage = zeros(1,length(subfiles));

%Number of artifact channels
artifact_chan = zeros(1,length(subfiles));

%proficiency analysis
prof_thresh = [70.0,70.0,40.0];
%found_prof: subs by pdm
found_prof = zeros(metadata.subjects,3);
%sub_prof: run where proficiency occurs-100 if never happens >66
sub_prof = 100.*ones(metadata.subjects,3);


%%

%loop through subjects
for sub_num = 1:metadata.subjects
    for sess_num = 1:11
        try
            file_fn = ['S',num2str(sub_num),'_Session_',num2str(sess_num),'.mat'];
            fprintf('loading file %s...\n',file_fn)
            tmp = load([path.data,file_fn]);
            BCI = tmp.BCI;
            
            file = file + 1;
            
            %get position information
            trialdata = BCI.TrialData;
            %sub_num = trialdata(1).subjectnumber;
            %sess_num = trialdata(1).sessionnumber;
            
            %add easy info to meta data
            runs = unique([trialdata.runnumber]);
            metadata.runs = metadata.runs + length(runs);
            metadata.trials = metadata.trials + length(trialdata);
            metadata.feedback_time = metadata.feedback_time + sum([trialdata.triallength]);
            sess_time = cellfun(@(x) size(x,2), BCI.time);
            metadata.total_time = metadata.total_time + sum(sess_time./1000);
            
            %add session artifact percentage
            artifact_percentage(file) = nanmean([trialdata.artifact]);
            
            %add number of artifact channels
            artifact_chan(file) =length(BCI.chaninfo.noisechan);
            
            %get pvc by run for subject/sess/pdm
            run_pvc = NaN(6,3);
            bci_table = struct2table(trialdata);
            for pdm = 1:3
                pdm_ind = find(bci_table.tasknumber == pdm);
                pdm_runs = unique(bci_table{pdm_ind,'runnumber'});
                for run = 1:length(pdm_runs)
                    run_ind = find(bci_table.runnumber == pdm_runs(run));
                    run_table = bci_table(run_ind,:);
                    run_pvc(run,pdm) = nanmean(run_table.result).*100;
                end
            end
            
            %add session to pvc mat
            avg_pvc = nanmean(run_pvc,1);
            bci_pvc(sub_num,sess_num,:) = avg_pvc;
            
            %check for proficiency
            for pdm = 1:3
                if ~found_prof(sub_num,pdm)
                    %3 checks, first block second block avg
                    if nanmean(run_pvc(1:3,pdm)) > prof_thresh(pdm)
                        found_prof(sub_num,pdm) = 1;
                        sub_prof(sub_num,pdm) = (sess_num - 1)*6 + 3;
                    elseif nanmean(run_pvc(4:6,pdm)) > prof_thresh(pdm)
                        found_prof(sub_num,pdm) = 1;
                        sub_prof(sub_num,pdm) = (sess_num - 1)*6 + 6;
                    elseif nanmean(run_pvc(:,pdm)) > prof_thresh(pdm)
                        found_prof(sub_num,pdm) = 1;
                        sub_prof(sub_num,pdm) = (sess_num - 1)*6 + 6;
                    end
                    
                end%if proficiency not found
            end
        catch
            warning('no file %s',file_fn)
        end%try to open
    end%sess num
end%sub num


%percent artifacts
figure
subplot(1,2,1)
hist(artifact_chan)
xlim([0,8])
xlabel('Number of noisy channels')
ylabel('Number of sessions')
subplot(1,2,2)
hist(artifact_percentage.*100,20)
xlim([0,8])
ylabel('Number of sessions')
xlabel('Percent of trials with artifacts(%)')
%Data analysis

%bci performance: classification accuracy
color_cell{1} = [0,0.45,0.74];
figure('position',[190 567 1050 520])
pdms = {'LR','UD','2D'};
for pdm = 1:3
    subplot(1,3,pdm)
    hold on
    plot(1:11,nanmean(bci_pvc(:,:,pdm),1),'Marker','.','MarkerSize',20,'MarkerFaceColor',color_cell{1},'LineStyle','none')
    hold on
    shadedErrorBar(1:11,nanmean(bci_pvc(:,:,pdm),1),nanstd(bci_pvc(:,:,pdm),[],1)./sqrt(size(bci_pvc,1)),{'markerfacecolor',color_cell{1},'color',color_cell{1}, 'LineWidth',1.5,'LineStyle','-'}, 1)
    xlim([1,11])
    xlabel('Session')
    ylabel('Percent valid correct (%)')
    %set(gca,'FontSize',20,'FontName','Myriad Pro')
    box off
    title(pdms{pdm})
end


%bci performance: proficiency plot
figure('position',[190 567 1050 520])
for pdm = 1:3
    [f,x] = ecdf(sub_prof(:,pdm)); %sub_prof = run when reach proficiency,f = % of
    subplot(1,3,pdm)
    stairs(x./6,f.*100)
    xlim([0,11.5])
    ylim([0,100])
    xlabel('Session')
    ylabel('Percent of proficient participants (%)')
end
%population, x = run when it happens

end%function