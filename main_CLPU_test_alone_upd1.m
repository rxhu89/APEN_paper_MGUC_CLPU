%% this is for testing the accuracy of the adaptive CLPU model alone, namely with the given on-off command
%% pre-otg is normal status: 
%% statu_ini, temp_ini @ norm; Uhvac_sched_ini = 1;  dp_ini = 0; khvac_ini = hvac_knorm_all(1)
%% here uses the optimization method(based on the adaptive CLPU formulation <MILP> for MGUC integration) to get the CLPU estimation (Yalmip+clpex), not the analytical method
%% add the case with constant otg but the time-varying temperature

clc;clear;close all

code_start_time = tic;
startclock = datestr(datetime('now'));

%% select simulation time
N_hour = 26;
step_time1 = 30;            % time step of the first stage, 30min

temprature_file_name ='real'; % select temperature REAL profile
% temprature_file_name = 'setup'; % selecte temperature (constant profile)

%% input folder and output folder
input_fd_hvac = strcat(pwd,'\','Input_Data_HVAC_pecan','\');
input_fd_hvac_sched = strcat(pwd,'\','Input_HVAC_sched','\');

output_data_fd_name = strcat( 'T_',temprature_file_name );
output_fd_hvac = strcat(pwd,'\Output_Data\',output_data_fd_name,'\'); % results
output_fd_hvac_est_step = strcat(pwd,'\Output_Data\',output_data_fd_name,'\HVAC_load_est_step\'); % the result of each steps est
output_fd_hvac_simu_step = strcat(pwd,'\Output_Data\',output_data_fd_name,'\HVAC_load_simu_step\'); % the hvac simulation

%% import temperature or set 
if strcmp(temprature_file_name, 'real')
    Tout_hour_all = readtable(strcat(input_fd_hvac,'Tout_from_mat_add','.xlsx'),'sheet', 'hour','ReadVariableNames',false); % [27;26;....]
else 
    Tout_hour_all = readtable(strcat(input_fd_hvac,'Tout_set_my','.xlsx'),'sheet', 'hour','ReadVariableNames',false); % [27;26;....]
end
Tout_hour_all = Tout_hour_all{1:N_hour,:};
N_Tout_profile = size(Tout_hour_all,2);
% Tout_hour = Tout_hour_all{1:N_hour,:};

%% interprate temperature data
% Tout_hour = Tout_hour{:,:};
    % 30_minute interporation
Tout_hour_interp = [Tout_hour_all;Tout_hour_all(end,:)] ;
interp_t1 = 1:length(Tout_hour_all)+1; % append a new data point to the end
interp_t2 = 1:1/2:length(Tout_hour_all)+1; % append a new data point to the end
Tout_30min0 = interp1(interp_t1,Tout_hour_interp,interp_t2 );
Tout_30min_all = 0*Tout_30min0(1:end-1,:);
for i = 1:length(Tout_30min0)-1
    Tout_30min_all(i,:) = (Tout_30min0(i,:)+Tout_30min0(i+1,:))/2 ; % the average
end  
% Tout_30min(Tout_30min>0) = 36; % mannually set temperature for tecting
    % 1_minute interporation
interp_t2 = 1:1/60:length(Tout_hour_all)+1; % append a new data point to the end
Tout_minute_all = interp1(interp_t1,Tout_hour_interp,interp_t2 );
Tout_minute_all = Tout_minute_all(1:end-1,:);

N_step_30min = length(Tout_30min_all);
N_1min = 30*N_step_30min;

%% set HVAC load on/off command Ug
hvac_sched_all =  readtable(strcat(input_fd_hvac_sched,'hvac_sched_senarios','.xlsx'),'ReadVariableNames',false); % [step,N_scene] 'sheet', 'hour',
hvac_sched_all = hvac_sched_all{:,:};
N_hvac_sched = size(hvac_sched_all,2); % with 30-min interval
% Uhvac_sched = hvac_sched_all(:,1);

% const outage duration cases
hvac_sched_all_const_otg =  readtable(strcat(input_fd_hvac_sched,'hvac_sched_senarios_const_otg_2h','.xlsx'),'ReadVariableNames',false); % [step,N_scene] 'sheet', 'hour',
hvac_sched_all_const_otg = hvac_sched_all_const_otg{:,:};
% Uhvac_sched = hvac_sched_all(:,1);

%% import adaptive HVAC data
% --- stage 1, based on group-ph
hvac_factor_base = readtable(strcat(input_fd_hvac,'factor_summary','.xlsx'), 'sheet', 'summary','PreserveVariableNames',1); % [G1-a;G1-b;G1-c;G2-a;G2-b....]

Tout_test = hvac_factor_base{:,1};
k_peak_dura0 = hvac_factor_base{:,2}; % delta_hour-Tout
k_decay_dura0 = hvac_factor_base{:,3}; % pu/hour
hvac_Pnorm_all0 = hvac_factor_base{:,4}; % kW
hvac_Peak_all0 = hvac_factor_base{1,5}; % kW
Dpeak_dura0 = hvac_factor_base{:,6}; % h

%% import hvac simulation para (RCQ)
N_hvac = 1100;
hvac_para_random = importdata(strcat(input_fd_hvac,'hvac_para_pecan_1100house_list','.mat')); % Tset ~24, TDB =1, diversified initial values, for the beginnital of the whole simulation
% Phvac_rate = hvac_para_random.P_rate(1:N_hvac);   % [1,N_house]  % power of HVAC, W
% hvac_para_mark = hvac_para_random.para_mark;  % 1 reasonable parameters
% Phvac_rate_total = Phvac_rate.*hvac_para_mark; % in line with hvac_Peak_all0
% Phvac_rate_total = sum(Phvac_rate_total)/1000;
temp_ini0 = hvac_para_random.x_ini(:,1:N_hvac);
stat_ini0 = hvac_para_random.stat_ini(1:N_hvac);

%% ****************** loop for all cases
for i = 1:N_Tout_profile
    Tout_hour = Tout_hour_all(:,i);
    Tout_30min = Tout_30min_all(:,i);
    Tout_minute = Tout_minute_all(:,i);
    
    %% get the adptive CLPU parameters according to the temperature
    k_peak_dura = zeros(N_step_30min,1); % same for all group&ph, under each Tout
    k_decay_dura = zeros(N_step_30min,1);% same for all group&ph, under each Tout
    hvac_knorm_all = zeros(N_step_30min,1);  % same for all group&ph, under each Tout
    Dpeak_dura = zeros(N_step_30min,1); % same for all group&ph, under each Tout
    for t =1:length(Tout_30min)
        Tout_value = round( Tout_30min(t) ); % only get the parameter table with 1 c interval
        Tout_ind = find( Tout_test == Tout_value);
        k_peak_dura(t) = k_peak_dura0(Tout_ind); % delta_hour-Tout
        k_decay_dura(t) = k_decay_dura0(Tout_ind); % pu/hour
        hvac_knorm_all(t) = hvac_Pnorm_all0(Tout_ind)/hvac_Peak_all0;
        Dpeak_dura(t) = Dpeak_dura0(Tout_ind); % x hour, per unit, 0.5step(hour) -> 30min
    end
    hvac_kpeak_dura_all = k_peak_dura;     % 30-min per step !!!!!!!!!!!!! SP
    hvac_kdecay_dura_all = k_decay_dura/2; % 30-min per step
    hvac_dpeak_dura_all = Dpeak_dura*2;    % 1hour to each step(30min)
    %%
    for j = 1 : N_hvac_sched
        Uhvac_sched = hvac_sched_all(:,j);
        Uhvac_sched_const_otg = hvac_sched_all_const_otg(:,j);
        if strcmp(temprature_file_name, 'real')
            mark_case = strcat('R',num2str(i),'_','U',num2str(j));
            mark_case_norm = strcat('R',num2str(i));
            mark_case_title = strcat('R',num2str(i),'-','U',num2str(j));
            mark_case_norm_title = strcat('R',num2str(i));
        else % using mannually temperature profile setup; constant or @ not real
            mark_case = strcat('S',num2str(i),'_','U',num2str(j));
            mark_case_norm = strcat('S',num2str(i));
            mark_case_title = strcat('S',num2str(i),'-','U',num2str(j));
            mark_case_norm_title = strcat('S',num2str(i));
        end
        
        fprintf(['\n','@@@@ Processing Case: ',mark_case, '\n',])

        %% adaptive CLPU optimization
        % initial value for CLPU
        Uhvac_sched_ini = 1;     % the initial status
%         dp_ini    = 0.3;   % peak increament_1 since initial off for an hour
%         dp_ini = hvac_dpeak_dura_all(1); % assume has lost diversity totally
        dp_ini = 0; % 
        dpeak_ini = dp_ini;             % peak increament_2
        dre_ini   = 0;    % the Peak_
        khvac_ini = hvac_knorm_all(1);  % clpu add part factor
        M_big = 50;
        %% ***************** modelling **************************************
        
        % == define all variables =
        yalmip('clear')  % https://github.com/yalmip/YALMIP/discussions/956; yalmip won't clear the previous same-name variables in loop like matlab
        % sdpvar(n,n) % SYMMETRIC!
        x_d       = sdpvar(N_step_30min,1,'full');  % d_t
        x_dpeak   = sdpvar(N_step_30min,1,'full');  % dpeak_t
        x_Usatu   = binvar(N_step_30min,1,'full'); % Usatu_t
        x_dre     = sdpvar(N_step_30min,1,'full');  % dre_t
        x_Udecay  = binvar(N_step_30min,1,'full'); % Udecay_t
        x_khvac   = sdpvar(N_step_30min,1,'full');   % khvac
        x_kclpu   = sdpvar(N_step_30min,1,'full');  % kclpu
        x_Pclpu   = sdpvar(N_step_30min,1,'full');  % kclpu
        % ==== contstraints =====
        constraints = [];
        % peak
        constraints = [constraints, ( 0 <= x_d <= M_big*(1-Uhvac_sched) ):'Peak_1'];
        constraints = [constraints, ( x_d(1) + M_big*Uhvac_sched(1) >= dp_ini + hvac_kpeak_dura_all(1)*(1-Uhvac_sched(1)) ):'Peak_2ini' ];
        if N_step_30min >1
            constraints = [constraints, ( x_d(2:end) + M_big*Uhvac_sched(2:end)  >= x_d(1:end-1) + hvac_kpeak_dura_all(2:end).*(1-Uhvac_sched(2:end)) ):'Peak_2' ];
        end
        % peak saturation
        constraints = [constraints, ( x_dpeak >= x_d - M_big*x_Usatu ):'Peak_satu1'];
        constraints = [constraints, ( x_dpeak >= hvac_dpeak_dura_all.*x_Usatu ):'Peak_satu2'];
        constraints = [constraints, ( x_Usatu <= 1 - Uhvac_sched ):'Peak_satu3'];
        % remaining peak
        constraints = [constraints, ( x_dre(1) - dre_ini >= dpeak_ini - Uhvac_sched_ini - M_big*(1-Uhvac_sched(1)) ):'Peak_re1' ];
        if N_step_30min >1
            constraints = [constraints, ( x_dre(2:end) - x_dre(1:end-1) >= x_dpeak(1:end-1) - Uhvac_sched(1:end-1) - M_big*(1-Uhvac_sched(2:end)) ):'Peak_re2' ];
        end
        constraints = [constraints, ( 0 <= x_dre <= M_big*Uhvac_sched ):'Peak_re3' ];
        % decay
        constraints = [constraints, ( -M_big*x_Udecay + 0.001*Uhvac_sched <= x_dre ): 'decay_1'];
        constraints = [constraints, ( x_dre <= M_big*(1-x_Udecay) ): 'decay_2'];
        constraints = [constraints, ( x_Udecay <= Uhvac_sched ): 'decay_3'];
        % hvac factor
        constraints = [constraints, ( x_khvac(1) >= khvac_ini - hvac_kdecay_dura_all(1)*x_Udecay(1) + (Uhvac_sched(1) - Uhvac_sched_ini) ): 'khvac_1'];
        if N_step_30min >1
            constraints = [constraints, ( x_khvac(2:end) >= x_khvac(1:end-1) - hvac_kdecay_dura_all(2:end).*x_Udecay(2:end) + Uhvac_sched(2:end) - Uhvac_sched(1:end-1) ): 'khvac_2'];
        end
        constraints = [constraints, (  hvac_knorm_all.*Uhvac_sched <= x_khvac <= M_big*Uhvac_sched ): 'khvac_3'];
        % kclpu factor
        constraints = [constraints, (  x_kclpu == x_khvac - hvac_knorm_all.*Uhvac_sched ): 'kclpu_1'];
        constraints = [constraints, (  x_Pclpu == hvac_Peak_all0 * x_khvac ): 'Pclpu_1'];
        % bounds
        constraints = [constraints, ( 0 <= x_dpeak <= M_big ):'Bound_1'];
        constraints = [constraints, ( 0 <= x_dre <=   M_big ):'Bound_2'];
        constraints = [constraints, ( 0 <= x_kclpu <= 1 ):'Bound_3'];
        % ===== objective =====
        objective = sum(sum(x_Pclpu)) + sum(sum(x_dpeak)) + sum(sum(x_dre)); % minimize
        % =======solve =============================================
        % https://yalmip.github.io/Extracting   https://yalmip.github.io/command/sdpsettings/
        options = sdpsettings('solver','cplex','debug',1,'verbose',2,'savesolverinput',1,'savesolveroutput',1);  % verbose: display,0: no; 1 simple-display; 2 full-display

        % options.cplex.mip.tolerances % https://groups.google.com/g/yalmip/c/6sbthHn0SaQ
        solver_time_max_set = 30; % 1500
        options.cplex.timelimit = solver_time_max_set;

        % mip_gap = 0.0005; % 0-1 default 0.0001 relative tolerance on the gap between the best integer objective  
        % options.cplex.mip.tolerances.mipgap = mip_gap; 

        fprintf(['\n','******** Beginning solve ','\n',])

        sol = optimize([constraints],objective,options);

        fval = value(objective);
        fprintf(['\n','******** fval: ',num2str(fval),'(kWh)\n',])
        % ===== Analyze error flags =====
        mark_get_sol = 0;
        if sol.problem == 0 || sol.problem == 3 % 0 Successfully solved; 1 Infeasible problem; 2 Unbounded objective function; 3 Maximum iterations exceeded
            redu_sd = 0;
            fprintf(['******** fval: ',num2str(fval),'\n',])
            fprintf(['   ^^^^ ',sol.info,'\n']) 
            mark_get_sol = 1;
        else  % if NO solution, msd  need to be adjusted, continute to resolve the optimization
            disp('Hmm, something went wrong!');    
            fprintf(['\n',' ^^^^ ',sol.info,'\n'])      
            error('!!!! Failure !!!! ')
        end
        %% get the results
        sol_d       = value(x_d); 
        sol_dpeak   = value(x_dpeak);
        sol_dre     = value(x_dre);
        sol_Usatu   = value(x_Usatu);
        sol_Udecay  = value(x_Udecay); 
        sol_khvac = value(x_khvac);
        sol_khvac(sol_khvac>=0.999) = 1; % !!!!
        sol_kclpu = value(x_kclpu);
        sol_Pclpu = value(x_Pclpu);
        sol_Phvac = sol_khvac*hvac_Peak_all0;
        sol_info = sol;

        est_Phvac_otg_adapt.Uhvac_sched_ini = Uhvac_sched_ini;
        est_Phvac_otg_adapt.dp_ini = dp_ini;
        est_Phvac_otg_adapt.dpeak_ini  = dpeak_ini;
        est_Phvac_otg_adapt.dre_ini = dre_ini;
        est_Phvac_otg_adapt.khvac_ini = khvac_ini;
        est_Phvac_otg_adapt.M_big = M_big;

        est_Phvac_otg_adapt.sol_d = sol_d;
        est_Phvac_otg_adapt.sol_dpeak = sol_dpeak;
        est_Phvac_otg_adapt.sol_dre = sol_dre;
        est_Phvac_otg_adapt.sol_Usatu = sol_Usatu;
        est_Phvac_otg_adapt.sol_Udecay = sol_Udecay;
        est_Phvac_otg_adapt.sol_khvac = sol_khvac;
        est_Phvac_otg_adapt.sol_kclpu = sol_kclpu;
        est_Phvac_otg_adapt.sol_Pclpu = sol_Pclpu;
        est_Phvac_otg_adapt.sol_Phvac = sol_Phvac;
        est_Phvac_otg_adapt.sol_info = sol_info;

        % post-process @ get the linear interpol of 5-min HVAC load for the optimized adaptive CLPU-enhanced MGUC
        sol_Phvac_5min = zeros(N_step_30min*6,1);
        khvac_step_5min_ini = 0; % initial value for 5-min resol (previously)
        for s = 1:N_step_30min
            Tout_30min_step = Tout_30min(s);
            hvac_knorm_step = hvac_knorm_all(s);
            hvac_kdecay_dura_step = hvac_kdecay_dura_all(s);
            
            if s == 1
                khvac_2side_step = [khvac_ini;sol_khvac(1)];
                Uhvac_sched_2side_step = [Uhvac_sched_ini; Uhvac_sched(1)] ; % previous and current step
                dre_2side_step = [dre_ini; sol_dre(1)] ;
            else % only 1 step left
                khvac_2side_step = sol_khvac(s-1:s,1); % previous and current step
                Uhvac_sched_2side_step = Uhvac_sched(s-1:s,1);
                dre_2side_step = sol_dre(s-1:s,1); % previous and current step; @ at the biginning of the current step
            end
            Usatu_step = sol_Usatu(s,1);
            Udecay_step = sol_Udecay(s,1);
            [sol_khvac_step_5min] = sub_HVAC_estimation_to_5min(s, mark_case, hvac_knorm_step, hvac_kdecay_dura_step, Uhvac_sched_2side_step,...
                                                                     khvac_2side_step, dre_2side_step, Usatu_step, Udecay_step, khvac_step_5min_ini);
            sol_Phvac_5min((s-1)*6+1:s*6,1) = hvac_Peak_all0*sol_khvac_step_5min;
            % update the 
            khvac_step_5min_ini = sol_khvac_step_5min(end);
        end
        est_Phvac_otg_adapt.sol_Phvac_5min = sol_Phvac_5min;
        save( strcat(output_fd_hvac, 'est_Phvac_otg_adapt','_',mark_case,'.mat'), 'est_Phvac_otg_adapt')
            

        %% %%%%%%%%%%%%%%%%%% HVAC load simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %% initialize hvacs for normal operation and for CLPU simulation; simulate all periods together
        N_hour_ini = 60; % the hours for initialization
        Tout_set_min_ini = Tout_minute(1)*ones(N_hour_ini*60,1); % simulate all periods together
        N_step_ini = length(Tout_set_min_ini);

        % HVAC status for initialization
        s = 0; % mark the rolling simulaiton
        Uhvac_norm_ini = 1;  %% on
        % mark_case_ini = strcat('ini','_','norm');
        [Phvac_house_ini_norm, Troom_house_ini_norm, status_end_ini_norm, temp_end_ini_norm] ...
            = HVAC_load_group_simu(s, Tout_set_min_ini, N_hvac, stat_ini0, temp_ini0, Uhvac_norm_ini, input_fd_hvac,output_fd_hvac_simu_step); % mark_case,


        Uhvac_otg_ini  = 0; %% off
        % mark_case_ini = strcat('ini','_','otg');
        [Phvac_house_ini_otg, Troom_house_ini_otg, status_end_ini_otg, temp_end_ini_otg] ...
            = HVAC_load_group_simu(s, Tout_set_min_ini, N_hvac, stat_ini0, temp_ini0, Uhvac_otg_ini, input_fd_hvac, output_fd_hvac_simu_step); % mark_case,
        % save data
        simu_Phvac_ini.Tout_set_min_ini = Tout_set_min_ini;
        simu_Phvac_ini.Uhvac_norm_ini = Uhvac_norm_ini;
        simu_Phvac_ini.Phvac_house_ini_norm = Phvac_house_ini_norm;
        simu_Phvac_ini.Troom_house_ini_norm = Troom_house_ini_norm;
        simu_Phvac_ini.status_ini_norm = status_end_ini_norm;
        simu_Phvac_ini.temp_ini_norm = temp_end_ini_norm;

        simu_Phvac_ini.Uhvac_norm_otg = Uhvac_otg_ini;
        simu_Phvac_ini.Phvac_house_ini_otg = Phvac_house_ini_otg;
        simu_Phvac_ini.Troom_house_ini_otg = Troom_house_ini_otg;
        simu_Phvac_ini.status_ini_otg = status_end_ini_otg;
        simu_Phvac_ini.temp_ini_otg = temp_end_ini_otg;
        save( strcat(output_fd_hvac, 'simu_Phvac_ini','_',mark_case,'.mat'), 'simu_Phvac_ini')

        %% normal load simulation
        Uhvac_norm = 1;
        % mark_case_norm = strcat('norm');
        [Phvac_house_norm, Troom_house_norm, status_end_norm, temp_end_norm]...
            = HVAC_load_group_simu(s, Tout_minute, N_hvac, status_end_ini_norm, temp_end_ini_norm, Uhvac_norm, input_fd_hvac, output_fd_hvac_simu_step); % mark_case,
        %  save data
        simu_Phvac_norm.Tout_minute  = Tout_minute;
        simu_Phvac_norm.status_ini_norm  = status_end_ini_norm;
        simu_Phvac_norm.temp_ini_norm  = temp_end_ini_norm;
        simu_Phvac_norm.Phvac_house_norm  = Phvac_house_norm;
        simu_Phvac_norm.Troom_house_norm  = Troom_house_norm;
        save( strcat(output_fd_hvac, 'simu_Phvac_norm','_',mark_case,'.mat'), 'simu_Phvac_norm')

        %% outage load simulation
        % mark_case_otg_real = strcat('otg_real');
        Phvac_house_otg_real = zeros(N_hvac,N_1min);
        Troom_house_otg_real = zeros(N_hvac,N_1min);
         for s = 1:N_step_30min
             Tout_set_min = Tout_minute((s-1)*30+1:s*30);
             Uhvac_otg = Uhvac_sched(s);
             if s == 1
%                  status_ini_otg = status_end_ini_otg;
%                  temp_ini_otg   = temp_end_ini_otg;
                 status_ini_otg = status_end_ini_norm;
                 temp_ini_otg   = temp_end_ini_norm;
             else
                 status_ini_otg = status_end_otg;
                 temp_ini_otg   = temp_end_otg;   
             end
             [Phvac_house_otg_each, Troom_house_otg_each, status_end_otg, temp_end_otg]...
                = HVAC_load_group_simu(s, Tout_set_min, N_hvac, status_ini_otg, temp_ini_otg, Uhvac_otg, input_fd_hvac,output_fd_hvac_simu_step); % mark_case,
            Phvac_house_otg_real(:,(s-1)*30+1:s*30) = Phvac_house_otg_each;
            Troom_house_otg_real(:,(s-1)*30+1:s*30) = Troom_house_otg_each;
         end
        simu_Phvac_otg_real.Tout_minute  = Tout_minute;
        simu_Phvac_otg_real.Uhvac_sched  = Uhvac_sched;
%         simu_Phvac_otg_real.status_ini_otg  = status_end_ini_otg;
%         simu_Phvac_otg_real.temp_ini_otg  = temp_end_ini_otg;
        simu_Phvac_otg_real.status_ini_otg  = status_end_ini_norm;
        simu_Phvac_otg_real.temp_ini_otg  = temp_end_ini_norm;
        simu_Phvac_otg_real.Phvac_house_otg_real  = Phvac_house_otg_real;
        simu_Phvac_otg_real.Troom_house_otg_real  = Troom_house_otg_real;
        save( strcat(output_fd_hvac, 'simu_Phvac_otg_real','_',mark_case,'.mat'), 'simu_Phvac_otg_real')


        %% use average temperature for CLPU simulation 
        mark_case_otg_Tave = strcat('otg_ave');
        Phvac_house_otg_Tave = zeros(N_hvac,N_1min);
        Troom_house_otg_Tave = zeros(N_hvac,N_1min);
        Tout_minute_Tave = mean(Tout_minute)*ones(N_1min,1);
         for s = 1:N_step_30min
             Tout_set_min_Tave = Tout_minute_Tave((s-1)*30+1:s*30);
             Uhvac_otg = Uhvac_sched(s);
             if s == 1
%                  status_ini_otg = status_end_ini_otg;
%                  temp_ini_otg   = temp_end_ini_otg;
                 status_ini_otg = status_end_ini_norm;
                 temp_ini_otg   = temp_end_ini_norm;
             else
                 status_ini_otg = status_end_otg;
                 temp_ini_otg   = temp_end_otg;   
             end
             [Phvac_house_otg_each, Troom_house_otg_each, status_end_otg, temp_end_otg]...
                = HVAC_load_group_simu(s, Tout_set_min_Tave, N_hvac, status_ini_otg, temp_ini_otg, Uhvac_otg, input_fd_hvac,output_fd_hvac_simu_step); % mark_case,
            Phvac_house_otg_Tave(:,(s-1)*30+1:s*30) = Phvac_house_otg_each;
            Troom_house_otg_Tave(:,(s-1)*30+1:s*30) = Troom_house_otg_each;
         end
        simu_Phvac_otg_Tave.Tout_minute_Tave  = Tout_minute_Tave;
        simu_Phvac_otg_Tave.Uhvac_sched     = Uhvac_sched;
%         simu_Phvac_otg_Tave.status_ini_otg  = status_end_ini_otg;
%         simu_Phvac_otg_Tave.temp_ini_otg  = temp_end_ini_otg;
        simu_Phvac_otg_Tave.status_ini_otg  = status_end_ini_norm;
        simu_Phvac_otg_Tave.temp_ini_otg  = temp_end_ini_norm;
        simu_Phvac_otg_Tave.Phvac_house_otg_Tave  = Phvac_house_otg_Tave;
        simu_Phvac_otg_Tave.Troom_house_otg_Tave  = Troom_house_otg_Tave;
        save( strcat(output_fd_hvac, 'simu_Phvac_otg_Tave','_',mark_case,'.mat'), 'simu_Phvac_otg_Tave')
        %-- special @ average temperature max (outage + decay) simulation

        %% use constant outage duration for CLPU simulation 
        mark_case_otg_Dconst = strcat('otg_ave_dura');
        Phvac_house_otg_Dconst = zeros(N_hvac,N_1min);
        Troom_house_otg_Dconst = zeros(N_hvac,N_1min);
         for s = 1:N_step_30min
             Tout_set_min = Tout_minute((s-1)*30+1:s*30);
             Uhvac_otg = Uhvac_sched_const_otg(s);
             if s == 1
%                  status_ini_otg = status_end_ini_otg;
%                  temp_ini_otg   = temp_end_ini_otg;
                 status_ini_otg = status_end_ini_norm;
                 temp_ini_otg   = temp_end_ini_norm;
             else
                 status_ini_otg = status_end_otg;
                 temp_ini_otg   = temp_end_otg;   
             end
             [Phvac_house_otg_each, Troom_house_otg_each, status_end_otg, temp_end_otg]...
                = HVAC_load_group_simu(s, Tout_set_min, N_hvac, status_ini_otg, temp_ini_otg, Uhvac_otg, input_fd_hvac,output_fd_hvac_simu_step); % mark_case,
            Phvac_house_otg_Dconst(:,(s-1)*30+1:s*30) = Phvac_house_otg_each;
            Troom_house_otg_Dconst(:,(s-1)*30+1:s*30) = Troom_house_otg_each;
         end
        simu_Phvac_otg_Tave.Tout_minute  = Tout_minute;
        simu_Phvac_otg_Tave.Uhvac_sched_const_otg     = Uhvac_sched_const_otg;
%         simu_Phvac_otg_Tave.status_ini_otg  = status_end_ini_otg;
%         simu_Phvac_otg_Tave.temp_ini_otg  = temp_end_ini_otg;
        simu_Phvac_otg_Dconst.status_ini_otg  = status_end_ini_norm;
        simu_Phvac_otg_Dconst.temp_ini_otg  = temp_end_ini_norm;
        simu_Phvac_otg_Dconst.Phvac_house_otg_Dconst  = Phvac_house_otg_Dconst;
        simu_Phvac_otg_Dconst.Troom_house_otg_Dconst  = Troom_house_otg_Dconst;
        save( strcat(output_fd_hvac, 'simu_Phvac_otg_Dconst','_',mark_case,'.mat'), 'simu_Phvac_otg_Dconst')
        %-- special @ average temperature max (outage + decay) simulation

        %% for adaptive CLPU @ HVAC load estimation
        simu_Phvac_norm_sum = sum(Phvac_house_norm);
        simu_Phvac_otg_real_sum = sum(Phvac_house_otg_real);
        simu_Phvac_otg_Tave_sum = sum(Phvac_house_otg_Tave);
        sol_Phvac_otg_est_adapt = kron(sol_Phvac,ones(30,1)); % to 1-min data
        sol_Phvac_otg_est_adapt_linear = kron(sol_Phvac_5min,ones(5,1)); % to 1-min data
        
        simu_Phvac_otg_Dconst_sum = sum(Phvac_house_otg_Dconst);

        %% plot HVAC & CLPU
%         fig_size   = [0.3,0.3,0.3,0.6]; % one day data
        fig_size   = [0.3,0.3,0.25,0.2]; % one day data
        xTick_time = [0:4:N_step_30min]./2;
        fg_1 = figure;
        set(gcf,'unit','normalized','position',fig_size)
%         sgtitle(['HVAC load comparison',' @ ', mark_case_title])
%         subplot(2,1,1)
        yyaxis left
%         area([1:N_1min]'/30/2,sol_Phvac_otg_est_adapt','FaceColor',[0.70,0.70,0.70],'facealpha',0.7,'EdgeColor','none','LineStyle','none') % Load 0.8 means 20% transparentl expected
%         hold on
        area([1:N_1min]'/30/2,sol_Phvac_otg_est_adapt_linear','FaceColor',[0.70,0.70,0.70],'facealpha',1,'EdgeColor',[0.50,0.50,0.50],'LineStyle','-','LineWidth',1) % 'none') % Load 0.8 means 20% transparentl expected
        hold on
        plot([1:N_1min]'/30/2,simu_Phvac_otg_real_sum,'r-','LineWidth',1);
        hold on
        plot([1:N_1min]'/30/2,simu_Phvac_otg_Tave_sum,'b--','LineWidth',1);
        hold on
        plot([1:N_1min]'/30/2,simu_Phvac_otg_Dconst_sum,'g--','LineWidth',1); % 'color',[0.39,0.83,0.07],
        hold on
        plot([1:N_1min]'/30/2,simu_Phvac_norm_sum,'k:','LineWidth',1.5);
        ylabel('HVAC Load (kW)')
        ylim([0,3000])
        set(gca,'YColor','k');
        
        yyaxis right
        plot([1:N_1min]'/30/2,Tout_minute,'-','LineWidth',1.5,'color',[1,0.7,0.0]) % [0.39,0.83,0.07]
        ylabel(strcat('Temperature (', char(0176),'C',')') )
        ylim([20,40])

        hLegend = legend({'Adaptive','Real','Fixed(temp)','Fixed(otg)','Normal','Temperature'},'Location','northwest','Orientation','horizontal'); %,'NumColumns',3); % vertical 'horizontal'
        set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent
        
        title(['HVAC load comparison',' @ ', mark_case_title])
%         xlabel('Hour')
        xlim([0,N_step_30min/2])    
        set(gca,'YColor','k');
        set(gca,'xTick',xTick_time) % [0:4:T_otg]./2)
        
%         subplot(2,1,2)
%         plot([1:N_1min]'/30/2,Troom_house_otg_real,'HandleVisibility','off')
%         hold on
%         plot([1:N_1min]'/30/2,Tout_minute,'-','LineWidth',1.5,'color',[1,0.7,0.0]) % [0.39,0.83,0.07]
%         hLegend = legend('Outdoor temperature','Location','northwest');
%         set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent
%         
%         title(['Real-otg Room temperature'])
%         ylabel(strcat('Temperature (', char(0176),'C',')') )
%         xlabel('Hour')
%         xlim([0,N_step_30min/2])
%         ylim([0,40])
%         set(gca,'xTick',xTick_time) % [0:4:T_otg]./2)
        
%         set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14); %,'FontWeight','Bold', 'LineWidth', 2);
        set(gca,'FontName','Times New Roman','FontSize',14)%,'LineWidth',1) set the axis & title & legend font
        set(gcf,'Color',[1,1,1])

        saveas(fg_1,strcat(output_fd_hvac, 'Phvac_comp','_',mark_case,'.fig'));        
        
    end 
    %% plot Hnormal HVAC power and temperature
    fig_size   = [0.3,0.3,0.3,0.6]; % one day data
    xTick_time = [0:4:N_step_30min]./2;
    fg_1 = figure;
    set(gcf,'unit','normalized','position',fig_size)
%         sgtitle(['HVAC load comparison',' @ ', mark_case_title])
    subplot(2,1,1)
    yyaxis left
    plot([1:N_1min]'/30/2,simu_Phvac_norm_sum,'k-','LineWidth',1.5);
    ylabel('HVAC Load (kW)')
    ylim([0,3000])
    set(gca,'YColor','k');

    yyaxis right
    plot([1:N_1min]'/30/2,Tout_minute,'-','LineWidth',1.5,'color',[1,0.7,0.0]) % [0.39,0.83,0.07]
    ylabel(strcat('Temperature (', char(0176),'C',')') )
    ylim([20,40])

    hLegend = legend('normal','Tout','Location','northwest'); % ,'Orientation','horizontal'); % vertical
    set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent

    title(['HVAC load normal',' @ ', mark_case_norm_title])
    xlabel('Hour')
    xlim([0,N_step_30min/2])    
    set(gca,'YColor','k');
    set(gca,'xTick',xTick_time) % [0:4:T_otg]./2)

    subplot(2,1,2)
    plot([1:N_1min]'/30/2,Troom_house_norm,'HandleVisibility','off')
    hold on
    plot([1:N_1min]'/30/2,Tout_minute,'-','LineWidth',1.5,'color',[1,0.7,0.0]) % [0.39,0.83,0.07]
    hLegend = legend('Outdoor temperature','Location','northwest');
    set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.6])) % [.5,.5,.5] is light gray; 0.8 means 20% transparent

    title(['Normal Room temperature'])
    ylabel(strcat('Temperature (', char(0176),'C',')') )
    xlabel('Hour')
    xlim([0,N_step_30min/2])
    ylim([0,40])
    set(gca,'xTick',xTick_time) % [0:4:T_otg]./2)

    set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14); %,'FontWeight','Bold', 'LineWidth', 2);
%         set(gca,'FontName','Times New Roman','FontSize',14)%,'LineWidth',1) set the axis & title & legend font
    set(gcf,'Color',[1,1,1])

    saveas(fg_1,strcat(output_fd_hvac, 'Phvac_norm','_',mark_case_norm,'.fig'));   
        
end   

%%

code_end_time = toc(code_start_time );
fprintf(['\n','!!!!!!!!! @ Finished All In: ',num2str(code_end_time),' sec.', '\n',])
