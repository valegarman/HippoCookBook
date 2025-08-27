%% 1- Automatically load all fiber and ripple data from all desired sessions and join them into a single file

% Directorio raíz (carpeta) donde están todas las subcarpetas por día (sessions) del animal a analizar
root_dir = 'Y:\unindexedSubjects\astro6';  %change this for each animal

% Listar todas las subcarpetas (ignorando '.' y '..')
folders = dir(root_dir);
folders = folders([folders.isdir]); 
folders = folders(~ismember({folders.name}, {'.', '..'}));

folder_names = {folders.name};

% Descartar las primeras 4 sessions (de prueba) y empezar con las sessions con behaviour. Customize this as desired for each animal
% folders = folders(5:end); 

% Preallocate cell arrays for joint data
all_green = {};     
all_ripples={};
all_iso={} ; 
all_red={} ; 
all_timestamps={} ; 
all_rec_mask={} ; 
all_rec_mask_timestamps={}; 
session_labels = {};  % para guardar el nombre (y orden) de las sessions


for i = 1:length(folders)  % Recorremos cada subcarpeta (session) de la carpeta global del animal:
    folder_path = fullfile(root_dir, folders(i).name);
    
    % Buscar archivo .mat deseado
    files_fiber = dir(fullfile(folder_path, '*FiberPhotometry.mat')); %buscar fiber data en esa session
    files_ripple= dir(fullfile(folder_path, '*ripples.events.mat')); %buscar ripple data en esa session
    
    if ~isempty(files_fiber) & ~isempty(files_ripple) %si se ha encontrado tanto fiber como ripple data:
        % Cargar files de interes
        fiber_data = load(fullfile(folder_path, files_fiber(1).name)); 
        ripple_data= load(fullfile(folder_path, files_ripple(1).name)); 
        
        % Extract the data of interest from those files
        green = fiber_data.fiber.green_PP.green_dFF_Smoothed;  
        iso= fiber_data.fiber.iso_PP.iso_dFF_Smoothed; %change to Smoothed_Z too when implemented in the other functions
        red= fiber_data.fiber.red_PP.red_dFF_Smoothed;
        timestamps=fiber_data.fiber.timestamps;
        rec_mask=fiber_data.fiber.events.maskSessions;
        rec_mask_timestamps=fiber_data.fiber.events.subSessions;
        ripples=ripple_data.ripples.peaks;
        
        % Guardar en el cell array que agrupará el data de todas las sessiones a medida que el loop itera
        all_green{end+1} = green(:);    
        all_iso{end+1} = iso(:);
        all_red{end+1} = red(:);
        all_timestamps{end+1} = timestamps(:);
        all_rec_mask{end+1} = rec_mask(:);
        all_rec_mask_timestamps{end+1} = rec_mask_timestamps(:,:);
        all_ripples{end+1} = ripples(:);
        session_labels{end+1} = folders(i).name; 
    else
        warning('FiberPhotometry.mat and/or ripples.mat not found in: %s', folder_path); %dar aviso si no se ha encontrado o fiber o ripple file en esa session
    end
end

% Opcional: mostrar cuántas señales se cargaron
fprintf('Se concatenaron %d señales.\n', length(all_green));

% Generate and save final boss fiber file containing the fiber-ripple data for all days
RipplesFiber.Data.green=all_green;
RipplesFiber.Data.iso=all_iso;
RipplesFiber.Data.red=all_red;
RipplesFiber.Data.timestamps=all_timestamps;
RipplesFiber.Data.rec_mask=all_rec_mask;
RipplesFiber.Data.rec_mask_timestamps=all_rec_mask_timestamps;
RipplesFiber.Data.ripples_peaks=all_ripples;
RipplesFiber.Data.session_labels=session_labels;

save(fullfile(root_dir, 'RipplesFiber.mat'), 'RipplesFiber');
disp(['Archivo guardado en: ', fullfile(root_dir, 'RipplesFiber.mat')]);



%% 2- Find the ripples associated to fiber signal

win=30; %seconds of signal we will look at
Basal=2; %seconds before ripple peak that we will use as basal for computing zscore
sorting=1; %seconds of fiber signal after ripple peak to use for computing the mean fiber response to ripple to sort them by intensity

for i=1:length(RipplesFiber.Data.green) %each session

    Mask_n=unique(RipplesFiber.Data.rec_mask{i}); %number of recs in the session (generally 3)
    Mask_use=RipplesFiber.Data.rec_mask{i}; %index of all recordings of the current session
    Session=['sess', extractAfter(RipplesFiber.Data.session_labels{i},'sess')]; %name of the session we are working on in the iteration
    Ripples_session=RipplesFiber.Data.ripples_peaks{i}; %timestamps of the ripple peaks identified for this session (across all recordings, including the ones without fiber)   

    for j=1:length(Mask_n) %each recording of the current session:   

        Rec=['rec',num2str(Mask_n(j))]; %name of the recording we are working on in the iteration
        rec_use=Mask_use==j; %index for the current recording of current session        

        Green_use=RipplesFiber.Data.green{i}(rec_use); %green signal from current rec of current session
        Red_use=RipplesFiber.Data.red{i}(rec_use); %red signal from current rec of current session
        Iso_use=RipplesFiber.Data.iso{i}(rec_use); %iso signal from current rec of current session
        Timestamps_use=RipplesFiber.Data.timestamps{i}(rec_use); %timestamps from current rec of current session

        %Find the ripples that fall under this rec window:
        Current_rec=find(rec_use);  %ID of the current recording of current session 
        Initial_use=Timestamps_use(1); %first timestamps of the rec
        Final_use=Timestamps_use(end); %last timestamps of the rec
        Ripples_use=Ripples_session;
        Ripples_use(Ripples_session<Initial_use+win)=[]; %Discard ripples happening before the start of this rec
        Ripples_use(Ripples_use>Final_use-win)=[]; %Discard ripples happening after the end of this rec
        [~, Ripple_ID]=ismember(Ripples_use,Ripples_session);
      
        %Now I kept only the ripples falling within the timestamps of the fiber recording of the iteration (with a 5s buffer to allow for baseline window)

        spacing=mean(diff(Timestamps_use(1:100))); %time separation between points
        sr=round(1/spacing); %find automatically the sampling rate (should always be the same, but just in case)
        Basal_use=(win-Basal)*sr:win*sr; %idx of signal before ripple peak that we will use as basal for computing z-scores
        win_sz=win*sr; %idx of signal we will look at in the plots (=time window)
        x_t=linspace(-win, win, (2*win_sz) + 1); %time vector for the plot

        %Preallocate variables that will store the fiber signal associated to each ripple occuring during that recording
        Green_ripple=[];
        Red_ripple=[];
        Iso_ripple=[];
        fiber_ripple_ts=[];
        
        for k=1:size(Ripples_use,1) %for each ripple in the rec
            current_ripple_ts=Ripples_use(k); %when did it happen
            [aa, current_fiber_ts]=min(abs((Timestamps_use-current_ripple_ts))); %find the corresponding fiber timestamps of the ripple peak           
            current_ts_vector= Timestamps_use(current_fiber_ts-(win_sz):current_fiber_ts+(win_sz)); %compute the time window we are interested on around the ripple
            fiber_ripple_ts=[fiber_ripple_ts,current_ts_vector]; %store the time vector into joint variable

            current_vector_green=Green_use(current_fiber_ts-(win_sz):current_fiber_ts+(win_sz)); %find the green signal values during the time around the ripple
            Green_ripple=[Green_ripple, current_vector_green]; %store the signal into joint variable

            current_vector_iso=Iso_use(current_fiber_ts-(win_sz):current_fiber_ts+(win_sz)); %same for iso
            Iso_ripple=[Iso_ripple, current_vector_iso];
    
            current_vector_red=Red_use(current_fiber_ts-(win_sz):current_fiber_ts+(win_sz)); %same for red
            Red_ripple=[Red_ripple, current_vector_red];

        end

        %Now, find the mean of the signal value for each timestamp across ripples. This generates the mean fiber response to a ripple
        Green_ripple_mean=mean(Green_ripple,2);
        Iso_ripple_mean=mean(Iso_ripple,2);
        Red_ripple_mean=mean(Red_ripple,2);   
                
        %With those means calculated, compute now z-scores using a baseline previously defined of 2s before ripple peak
        Green_ripple_basal=mean(Green_ripple_mean(Basal_use,:));
        Green_ripple_std=std(Green_ripple_mean(Basal_use,:));
        Green_ripple_zs=(Green_ripple_mean-Green_ripple_basal)./Green_ripple_std; %compute the z-score having previously done the mean fiber response to the ripple
                
        Iso_ripple_basal=mean(Iso_ripple_mean(Basal_use,:));
        Iso_ripple_std=std(Iso_ripple_mean(Basal_use,:));
        Iso_ripple_zs=(Iso_ripple_mean-Iso_ripple_basal)./Iso_ripple_std;
        
        Red_ripple_basal=mean(Red_ripple_mean(Basal_use,:));
        Red_ripple_std=std(Red_ripple_mean(Basal_use,:));
        Red_ripple_zs=(Red_ripple_mean-Red_ripple_basal)./Red_ripple_std;

                
        %Store the variables into the 'final boss' structure for later recall during stats and plotting
            %Store them by sessions     
        if j==1 %preallocate
            RipplesFiber.Analysis.Timestamps_Ripple.(Session)=[];
            RipplesFiber.Analysis.Green_Ripple.(Session).All_ripples=[];
            RipplesFiber.Analysis.Iso_Ripple.(Session).All_ripples=[];
            RipplesFiber.Analysis.Red_Ripple.(Session).All_ripples=[];
            RipplesFiber.Analysis.Green_Ripple.(Session).Rec_mask=[]; % To keep track of rec identity of each ripple of the session
            RipplesFiber.Analysis.Green_Ripple.(Session).Rec_mask=[RipplesFiber.Analysis.Green_Ripple.(Session).Rec_mask ; ones(size(Green_ripple,2),1)]; %mask for first rec of the session
            RipplesFiber.Analysis.Green_Ripple.(Session).RippleID=[]; % To keep track of ripple identity. With the rec mask, we know to which recording do these IDs belong to.             
        end

        RipplesFiber.Analysis.Green_Ripple.(Session).RippleID=[RipplesFiber.Analysis.Green_Ripple.(Session).RippleID ; Ripple_ID ];
        RipplesFiber.Analysis.Timestamps_Ripple.(Session)=[RipplesFiber.Analysis.Timestamps_Ripple.(Session), fiber_ripple_ts]; %stores the timestamps vector of fiber signal for each ripple across recs of the given session. (Column=ripple)
        
        RipplesFiber.Analysis.Green_Ripple.(Session).All_ripples=[RipplesFiber.Analysis.Green_Ripple.(Session).All_ripples, Green_ripple]; %Variable has, for the 601 timestamps (rows), the fiber values for all the ripples (columns) of that full session (so the 3 recs concatenated)
        RipplesFiber.Analysis.Iso_Ripple.(Session).All_ripples= [RipplesFiber.Analysis.Iso_Ripple.(Session).All_ripples, Iso_ripple]; 
        RipplesFiber.Analysis.Red_Ripple.(Session).All_ripples=[RipplesFiber.Analysis.Red_Ripple.(Session).All_ripples, Red_ripple]; 

        RipplesFiber.Analysis.Green_Ripple.(Session).Mean_zs(:,j)=Green_ripple_zs; %collapse ripple dimension => mean fiber response to a ripple in each of the 3 recs of that session (Column 1=HCF1, 2=Maze, 3=HCF2)
        RipplesFiber.Analysis.Iso_Ripple.(Session).Mean_zs(:,j)=Iso_ripple_zs;
        RipplesFiber.Analysis.Red_Ripple.(Session).Mean_zs(:,j)=Red_ripple_zs;

        if j==2
            RipplesFiber.Analysis.Green_Ripple.(Session).Rec_mask=[RipplesFiber.Analysis.Green_Ripple.(Session).Rec_mask ; ones(size(Green_ripple,2),1)+1]; %mask for second rec of the session
        end

        if j==3
             RipplesFiber.Analysis.Green_Ripple.(Session).Mean_zs_aR=mean(RipplesFiber.Analysis.Green_Ripple.(Session).Mean_zs,2); %collapse also rec dimension => mean fiber reponse to a ripple for that session. 'aR' = across recordings
             RipplesFiber.Analysis.Iso_Ripple.(Session).Mean_zs_aR=mean(RipplesFiber.Analysis.Iso_Ripple.(Session).Mean_zs,2);
             RipplesFiber.Analysis.Red_Ripple.(Session).Mean_zs_aR=mean(RipplesFiber.Analysis.Red_Ripple.(Session).Mean_zs,2);
             RipplesFiber.Analysis.Green_Ripple.(Session).Rec_mask=[RipplesFiber.Analysis.Green_Ripple.(Session).Rec_mask ; ones(size(Green_ripple,2),1)+2]; %mask for third rec of the session
        end        
        
        
            
            %Store them by recs        
        if i==1 %preallocate
            RipplesFiber.Analysis.Timestamps_Ripple.(Rec)=[];
            RipplesFiber.Analysis.Green_Ripple.(Rec).All_ripples=[];
            RipplesFiber.Analysis.Iso_Ripple.(Rec).All_ripples=[];
            RipplesFiber.Analysis.Red_Ripple.(Rec).All_ripples=[];   

            RipplesFiber.Analysis.Green_Ripple.(Rec).RippleID=[]; %these ripple IDs will not be consecutive since they belong to different days and therefore ID restarts. But with the session mask we can still trace their full identity since we know to what session do they belong and the ID within that session
            
            if j==1
                RipplesFiber.Analysis.n_Ripples=[];
            end

        end
        
        RipplesFiber.Analysis.Green_Ripple.(Rec).RippleID=[RipplesFiber.Analysis.Green_Ripple.(Rec).RippleID ; Ripple_ID ];
        RipplesFiber.Analysis.Timestamps_Ripple.(Rec)=[RipplesFiber.Analysis.Timestamps_Ripple.(Rec), fiber_ripple_ts]; %stores the timestamps vector of fiber signal for each ripple across sessions of the given rec. (Column=ripple)

        RipplesFiber.Analysis.Green_Ripple.(Rec).All_ripples=[RipplesFiber.Analysis.Green_Ripple.(Rec).All_ripples, Green_ripple]; %fiber signal during each ripple (column) across sessions for a given recording.
        RipplesFiber.Analysis.Iso_Ripple.(Rec).All_ripples=[RipplesFiber.Analysis.Iso_Ripple.(Rec).All_ripples, Iso_ripple ];
        RipplesFiber.Analysis.Red_Ripple.(Rec).All_ripples=[RipplesFiber.Analysis.Red_Ripple.(Rec).All_ripples, Red_ripple ];
        
        
        RipplesFiber.Analysis.Green_Ripple.(Rec).Mean_zs(:,i)=Green_ripple_zs; %collapse ripple dimension => mean fiber response to a ripple for a given recording in each session. Order of columns (=sessions) can be checked in RipplesFiber.Data.session_labels
        RipplesFiber.Analysis.Iso_Ripple.(Rec).Mean_zs(:,i)=Iso_ripple_zs;
        RipplesFiber.Analysis.Red_Ripple.(Rec).Mean_zs(:,i)=Red_ripple_zs;

        RipplesFiber.Analysis.n_Ripples=[RipplesFiber.Analysis.n_Ripples, size(Green_ripple,2)];

    end

end

%Reorganise the number ripples variable so that it matches the 3recs by 14 sessions structure
RipplesFiber.Analysis.n_Ripples=reshape(RipplesFiber.Analysis.n_Ripples, 3, length(RipplesFiber.Data.session_labels))'; %Now the numbers of ripples are indicated as a function of rec (column) and session (row)

%Now, compute the mean of each rec across sessions and add it to the plots
RipplesFiber.Analysis.Green_Ripple.rec1.Mean_zs_aS=mean(RipplesFiber.Analysis.Green_Ripple.rec1.Mean_zs,2); %Collapse session dimension => mean fiber signal during a ripple in rec1
RipplesFiber.Analysis.Green_Ripple.rec2.Mean_zs_aS=mean(RipplesFiber.Analysis.Green_Ripple.rec2.Mean_zs,2); %same for rec2
RipplesFiber.Analysis.Green_Ripple.rec3.Mean_zs_aS=mean(RipplesFiber.Analysis.Green_Ripple.rec3.Mean_zs,2); %same for rec2

RipplesFiber.Analysis.Iso_Ripple.rec1.Mean_zs_aS=mean(RipplesFiber.Analysis.Iso_Ripple.rec1.Mean_zs,2); %same for iso
RipplesFiber.Analysis.Iso_Ripple.rec2.Mean_zs_aS=mean(RipplesFiber.Analysis.Iso_Ripple.rec2.Mean_zs,2);
RipplesFiber.Analysis.Iso_Ripple.rec3.Mean_zs_aS=mean(RipplesFiber.Analysis.Iso_Ripple.rec3.Mean_zs,2);

RipplesFiber.Analysis.Red_Ripple.rec1.Mean_zs_aS=mean(RipplesFiber.Analysis.Red_Ripple.rec1.Mean_zs,2); %same for red
RipplesFiber.Analysis.Red_Ripple.rec2.Mean_zs_aS=mean(RipplesFiber.Analysis.Red_Ripple.rec2.Mean_zs,2);
RipplesFiber.Analysis.Red_Ripple.rec3.Mean_zs_aS=mean(RipplesFiber.Analysis.Red_Ripple.rec3.Mean_zs,2);

% Compute a mask to distinguish the ripples belonging to different session
RipplesFiber.Analysis.Green_Ripple.rec1.Session_mask=[];
RipplesFiber.Analysis.Green_Ripple.rec2.Session_mask=[];
RipplesFiber.Analysis.Green_Ripple.rec3.Session_mask=[];

for i=1:length(RipplesFiber.Data.session_labels)   
    RipplesFiber.Analysis.Green_Ripple.rec1.Session_mask=[RipplesFiber.Analysis.Green_Ripple.rec1.Session_mask ; zeros(RipplesFiber.Analysis.n_Ripples(i,1),1)+i];  
    RipplesFiber.Analysis.Green_Ripple.rec2.Session_mask=[RipplesFiber.Analysis.Green_Ripple.rec2.Session_mask ; zeros(RipplesFiber.Analysis.n_Ripples(i,2),1)+i];  
    RipplesFiber.Analysis.Green_Ripple.rec3.Session_mask=[RipplesFiber.Analysis.Green_Ripple.rec3.Session_mask ; zeros(RipplesFiber.Analysis.n_Ripples(i,3),1)+i];  
end

RipplesFiber.Analysis.plot_time=x_t;


%% 3 - FIGURES
%Note: in plot_ripple_align_alternative.m the code allows to plot the
%heatmaps either by ripple order (so time/session-wise) or by fiber
%response intensity to the ripple easily with a switch-case. In this code,
%to alternate between either option you have to comment and uncomment a
%couple lines
sorting_win=win*sr+1:(win+sorting)*sr+1;

for i=1:length(RipplesFiber.Data.session_labels)
    Session=['sess', extractAfter(RipplesFiber.Data.session_labels{i},'sess')];

    % 3.1: Traces and heatmaps for every rec of each session (just for green. If red and iso desired, run ripple_align.m or ripple_align_across_days.m or write the code changing 'green' for 'red' or 'iso')
    fig_rec(i)=figure('Name',[Session, ' - Recs']);  

    ax1=subplot(6,3,1);        
    plot(x_t,RipplesFiber.Analysis.Green_Ripple.(Session).Mean_zs(:,1),'color',[0.2 0.5 0.2]),hold on
    % xlim([-3 3])
    xlim([x_t(1) x_t(end)])
    title(['Rec1', ' - Mean trace'])
   
    ax2=subplot(6,3,[4,7,10,13,16]);
    idx_R1=RipplesFiber.Analysis.Green_Ripple.(Session).All_ripples(:,1:RipplesFiber.Analysis.n_Ripples(i,1));
    Green_single_Z_R1=(idx_R1-mean(idx_R1(Basal_use,:)))./std(idx_R1(Basal_use,:),[],1);
    Green_single_Z_R1_sort=mean(Green_single_Z_R1(sorting_win,:)); 
    % imagesc(x_t,1:size(Green_single_Z_R1, 2),Green_single_Z_R1'); %if by natural order
    imagesc_ranked(x_t,[],Green_single_Z_R1',[-3 3],Green_single_Z_R1_sort'); %if sorted by response intensity to ripple
    caxis([-3 3]);
    colormap(jet);          
    xlabel('Time');
    ylabel('Ripples');
    title(['Rec1', ' - Heatmap -  ', 'Num ripples = ' num2str(RipplesFiber.Analysis.n_Ripples(i,1))])
    pos1 = get(ax2, 'Position');
    pos1(2) = pos1(2) - 0.05; % sube un poco el primer subplot
    set(ax2, 'Position', pos1);


    ax3=subplot(6,3,2);
    plot(x_t,RipplesFiber.Analysis.Green_Ripple.(Session).Mean_zs(:,2),'color',[0.2 0.5 0.2]),hold on
    % xlim([-3 3])
    xlim([x_t(1) x_t(end)])
    title(['Rec2',' - Mean trace'])

    ax4=subplot(6,3,[5,8,11,14,17]);
    idx_R2=RipplesFiber.Analysis.Green_Ripple.(Session).All_ripples(:,RipplesFiber.Analysis.n_Ripples(i,1)+1:RipplesFiber.Analysis.n_Ripples(i,1)+RipplesFiber.Analysis.n_Ripples(i,2));
    Green_single_Z_R2=(idx_R2-mean(idx_R2(Basal_use,:)))./std(idx_R2(Basal_use,:),[],1); 
    Green_single_Z_R2_sort=mean(Green_single_Z_R2(sorting_win,:)); 
    % imagesc(x_t,1:size(Green_single_Z_R2, 2),Green_single_Z_R2'); 
    imagesc_ranked(x_t,[],Green_single_Z_R2',[-3 3],Green_single_Z_R2_sort'); %if sorted by response intensity to ripple
    caxis([-3 3]);
    colormap(jet);          
    xlabel('Time');
    ylabel('Ripples');
    title(['Rec2', ' - Heatmap -  ', 'Num ripples = ' num2str(RipplesFiber.Analysis.n_Ripples(i,2))])
    pos1 = get(ax4, 'Position');
    pos1(2) = pos1(2) - 0.05; % sube un poco el primer subplot
    set(ax4, 'Position', pos1);


    ax5=subplot(6,3,3);
    plot(x_t,RipplesFiber.Analysis.Green_Ripple.(Session).Mean_zs(:,3),'color',[0.2 0.5 0.2]),hold on
    % xlim([-3 3])
    xlim([x_t(1) x_t(end)])
    title(['Rec3', ' - Mean trace'])

    ax6=subplot(6,3,[6,9,12,15,18]);
    idx_R3=RipplesFiber.Analysis.Green_Ripple.(Session).All_ripples(:,RipplesFiber.Analysis.n_Ripples(i,1)+RipplesFiber.Analysis.n_Ripples(i,2)+1:end);
    Green_single_Z_R3=(idx_R3-mean(idx_R3(Basal_use,:)))./std(idx_R3(Basal_use,:),[],1);
    Green_single_Z_R3_sort=mean(Green_single_Z_R3(sorting_win,:)); 
    % imagesc(x_t,1:size(Green_single_Z_R3, 2),Green_single_Z_R3'); 
    imagesc_ranked(x_t,[],Green_single_Z_R3',[-3 3],Green_single_Z_R3_sort'); %if sorted by response intensity to ripple
    caxis([-3 3]);
    colormap(jet);       
    cb=colorbar(ax6);
    xlabel('Time');
    ylabel('Ripples');
    title(['Rec3', ' - Heatmap -  ', 'Num ripples = ' num2str(RipplesFiber.Analysis.n_Ripples(i,3))])
    
    pos1 = get(ax6, 'Position');
    pos1(2) = pos1(2) - 0.05; % sube un poco el primer subplot
    set(ax6, 'Position', pos1);
    pos_ax6 = get(ax6, 'Position');
    pos_cb = get(cb, 'Position');
    pos_cb(1) = pos_ax6(1) + pos_ax6(3) + 0.01;  % 0.01 de espacio extra a la derecha
    pos_cb(2) = pos_ax6(2);                        % misma altura que ax3g
    pos_cb(4) = pos_ax6(4);                        % misma altura que ax3g    
    set(cb, 'Position', pos_cb);
    
    sgtitle([Session, ' - Green'],'FontSize',16)



    % 3.2: Mean traces across recs for that session and heatmap of all ripples of the session 
    fig_rec_mean(i)=figure('Name',[Session, ' - Recs mean']);  

    ax7 = subplot(6,1,1);
    plot(x_t,RipplesFiber.Analysis.Green_Ripple.(Session).Mean_zs_aR,'k')
    xlim([x_t(1) x_t(end)])
    set(ax1, 'XTickLabel', []); 
    title(['Mean trace ', Session])

    ax8 = subplot(6,1,2:6);
    Green_single_Z_aR=[Green_single_Z_R1, Green_single_Z_R2, Green_single_Z_R3];
    Green_single_Z_aR_sort=mean(Green_single_Z_aR(sorting_win,:)); 
    % imagesc(x_t,1:size(Green_single_Z_aR, 2),Green_single_Z_aR'); 
    imagesc_ranked(x_t,[],Green_single_Z_aR',[-3 3],Green_single_Z_aR_sort'); %if sorted by response intensity to ripple
    caxis([-3 3]);
    colormap(jet);       
    colorbar;
    xlabel('Time');
    ylabel('Ripples');
    title(['Green - All Ripples in ' , Session])

    pos1 = get(ax7, 'Position');
    pos2 = get(ax8, 'Position');
    pos2(1) = pos1(1);         % igualar la posición x
    pos2(3) = pos1(3);         % igualar el ancho
    pos2(2) = pos2(2) - 0.05;
    set(ax8, 'Position', pos2);

    sgtitle(['Green - ' Session])

end



% 3.3: Heatmaps + individual and mean trace across sessions of each rec
sort_method=1;
for j=1:max(unique(RipplesFiber.Data.rec_mask{1, 1}))
    Rec=['rec',num2str(j)];
    Green_single_Z_R_aS=(RipplesFiber.Analysis.Green_Ripple.(Rec).All_ripples - mean(RipplesFiber.Analysis.Green_Ripple.(Rec).All_ripples(180:300,:)))./std(RipplesFiber.Analysis.Green_Ripple.(Rec).All_ripples(180:300,:),[],1); %do the z-score without previous mean

    switch sort_method
        case 1                 
            figure_aS(i)=figure('Name',['fig aS - ',Rec]); %sorted per fiber response to the 1st sec after ripple peak
            
            ax1 = subplot(6,1,1);
            hold on
            for k=1:length(RipplesFiber.Data.session_labels)
                plot(x_t,RipplesFiber.Analysis.Green_Ripple.(Rec).Mean_zs(:,k))
            end
            plot(x_t,RipplesFiber.Analysis.Green_Ripple.(Rec).Mean_zs_aS,'k','LineWidth',1.5)
            hold off
            xlim([x_t(1) x_t(end)])
            set(ax1, 'XTickLabel', []);    
            
            ax2 = subplot(6,1,2:6);
            Green_single_Z_R3_aS_sort=mean(Green_single_Z_R_aS(301:361,:)); 
            
            imagesc_ranked(x_t,[],Green_single_Z_R_aS',[-3 3],Green_single_Z_R3_aS_sort'); 
            colormap(jet); 
            colorbar;
            xlabel('Time');
            ylabel('Ripples');
            sgtitle(['Green - Ripples in ',Rec,' across all sessions'],'FontSize',12)
            
            pos1 = get(ax1, 'Position');
            pos2 = get(ax2, 'Position');
            pos2(1) = pos1(1);         % igualar la posición x
            pos2(3) = pos1(3);         % igualar el ancho
            set(ax2, 'Position', pos2);
            
        case 2    
            figure('Name','fig_hm_R3_aS') %sorted per session (so day-wise)

            ax1 = subplot(6,1,1);
            hold on
            for k=1:length(RipplesFiber.Data.session_labels)
                plot(x_t,RipplesFiber.Analysis.Green_Ripple.(Rec).Mean_zs(:,k))
            end
            plot(x_t,RipplesFiber.Analysis.Green_Ripple.(Rec).Mean_zs_aS,'k','LineWidth',1.5)
            hold off
            xlim([x_t(1) x_t(end)])
            set(ax1, 'XTickLabel', []); 
            
            ax2 = subplot(6,1,2:6);
            imagesc(x_t,1:size(Green_single_Z_R_aS, 2),Green_single_Z_R_aS'); 
            caxis([-3 3]);
            colormap(jet); 
            colorbar;
            
            %To indicate where does each session end in the heatmap:
            cumulativeRipples_rec = cumsum(RipplesFiber.Analysis.n_Ripples(:,j));
            for i = 1:length(RipplesFiber.Data.session_labels)-1
                yline(cumulativeRipples_rec(i), '--', 'LineWidth', 1.5);
            end
            
            xlabel('Time');
            ylabel('Ripples');
            sgtitle(['Green - Ripples in ',Rec,' across all sessions'],'FontSize',14)
            
            pos1 = get(ax1, 'Position');
            pos2 = get(ax2, 'Position');
            pos2(1) = pos1(1);         % igualar la posición x
            pos2(3) = pos1(3);         % igualar el ancho
            set(ax2, 'Position', pos2);
    end

end


%For plotting fiber trace around all ripples, by recs and mean

rec1=RipplesFiber.Analysis.Green_Ripple.sess1.All_ripples(:,RipplesFiber.Analysis.Green_Ripple.sess1.Rec_mask ==1);
figure()
for i=1:size(rec1,2)
    hold on; plot(x_t,rec1(:,i));
    hold off
end

rec2=RipplesFiber.Analysis.Green_Ripple.sess1.All_ripples(:,RipplesFiber.Analysis.Green_Ripple.sess1.Rec_mask ==2);
figure()
for i=1:size(rec2,2)
    hold on; plot(x_t,rec2(:,i));
    hold off
end

rec3=RipplesFiber.Analysis.Green_Ripple.sess1.All_ripples(:,RipplesFiber.Analysis.Green_Ripple.sess1.Rec_mask ==3);
figure()
for i=1:size(rec3,2)
    hold on; plot(x_t,rec3(:,i));
    hold off
end

figure()
hold on
plot(x_t,RipplesFiber.Analysis.Green_Ripple.sess1.Mean_zs  (:,1))
plot(x_t,RipplesFiber.Analysis.Green_Ripple.sess1.Mean_zs  (:,2))
plot(x_t,RipplesFiber.Analysis.Green_Ripple.sess1.Mean_zs  (:,3))
plot(x_t,RipplesFiber.Analysis.Green_Ripple.sess1.Mean_zs_aR,'k', 'LineWidth',1.5)
hold off



%% 4- Stats






%% 5 - Save results

save('RipplesFiber.mat', 'RipplesFiber')