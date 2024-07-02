% This Matlab code performs the molecular dynamics (MD) simulation for one-component Lennard-Jones (LJ) fluid. 
% A velocity-Verlet scheme is used in numerical integration. 
% Results for a potential energy per particle and virial pressure (beta*p/rho) from this calculation
% are compared with those obtained by Verlet, Ref. [2]. 
%
% This matlab code is written in similar way that is used in Fortran code in 
% the Appendix of Ref.[1]. 
%
% Ref. [1] D. Heermann, "Computer Simulation Methods in Theoretical Physics", 2nd edition, (1989);
% Ref. [2] L. Verlet, Phys. Rev. v159, p98 (1967); 
% Ref. [3] D. Frenkel and B. Smit, "Understanding Molecular Simulation", Acedmic Press (2002);
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
%
% July 2, 2024 & University of North Dakota 
%
function [] = classical_MD_for_LJ_fluid
clc; 
format short 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
number_of_dim = 4;           % numerical parameter & you may change it
npart = 4 * number_of_dim^3; % number of particle 
den = 0.8442;  %  reduced density 
T_ref = 0.728; %  reduced reference temperature
%
side = (npart/den)^(1/3);  % length of box side
A = side / number_of_dim;
%
sideh = side * 0.5;        
rc = 2.5;                  % cut-off parameter in the Lennard-Jones potential  
rc2 = rc * rc;
%
mass = 1.;  %  reduced mass 
dt = 0.010; %  reduced time-step 
%
[x, y, z] = coord_x_y_z(A, npart, number_of_dim);
%
%figure(1)
%scatter3(x,y,z, "blue", "filled")
[vx, vy, vz] = init_vel(npart, T_ref, mass);
%
[ax,ay,az] = force_acc(x, y, z, side, sideh, rc2, npart, mass);

%%%%%%%%%%%%%%%%%%%%%%%%
fileID_save_data_1 = fopen('classical_MD_for_LJ_fluid.txt','w');
%
Nstep = 10000.;
%
for ii = 1:Nstep
    %
    x = x + vx * dt + 0.5 * ax * dt * dt;  % velocity Verlet scheme: stage 1
    y = y + vy * dt + 0.5 * ay * dt * dt;    
    z = z + vz * dt + 0.5 * az * dt * dt;
    %
    [x, y, z] = periodic_BC(x, y, z, npart, side);
    %
    [ax1, ay1, az1, epot, vir] = force_acc(x, y, z, side, sideh, rc2, npart, mass); % velocity Verlet scheme: stage 2
    %
    vx = vx + 0.5 * (ax + ax1) * dt;
    vy = vy + 0.5 * (ay + ay1) * dt;
    vz = vz + 0.5 * (az + az1) * dt;
    %
    ax = ax1;
    ay = ay1;
    az = az1;
    %
    sm_v2 = sum(vx.*vx + vy.*vy + vz.*vz);
    %
    K = mass.*sm_v2/(2.*npart);
    T_inst = (2/3).*K;
    %
    % scaling factor , fs
    fs = sqrt(T_ref/T_inst);
    %
    vx = fs.*vx; 
    vy = fs.*vy; 
    vz = fs.*vz;

    %
    output = [ii*dt, epot, vir, T_inst];

    fprintf(fileID_save_data_1, '%4.4f \t %4.6f \t %4.6f \t %8.4f\n', output); 

end
fclose(fileID_save_data_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_md_data = fopen('classical_MD_for_LJ_fluid.txt', 'r');               % 
read_md_data = textscan(read_md_data, '%f %f %f %f');
md_step_ii = read_md_data{1};
md_epot = read_md_data{2};
md_vir = read_md_data{3};
md_ave_temp = read_md_data{4};


%%%
T_inst_ave = sum(md_ave_temp)/length(md_ave_temp); % instantaneous temperature computed by the MD simulation 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculation potential energy (U*/epsilon) per particle
%
Potential_energy_tail_correction_per_particle = ((8/3)*pi*den)*((1/3)*(1./rc^9) - (1./rc^3)); % from Ref. [3]. 
Potential_energy_per_partcile_without_tail_correction = sum(md_epot)/length(md_epot);
Potential_energy_per_partcile_with_tail_correction = sum(md_epot)/length(md_epot) + Potential_energy_tail_correction_per_particle;
%%%
[den, T_ref, Potential_energy_per_partcile_without_tail_correction, Potential_energy_per_partcile_with_tail_correction]

%[den, T_ref, Potential_energy_per_partcile_without_tail_correction, Potential_energy_per_partcile_with_tail_correction]
% N = 256
% 0.5426    1.4040   -3.3329   -3.6234 vs -3.63 from Ref. [2]
% 0.8500    2.2020   -4.2777   -4.7329 vs -4.76 from Ref. [2]
% 0.8800    0.9400   -5.5625   -6.0336 vs -5.84 from Ref. [2]
% 0.4500    1.7440   -2.6705   -2.9115 vs -2.90 from Ref. [2]
% 0.4000    1.4600   -2.4702   -2.6844 vs -2.72\pm 0.01 from Ref. [2]
%
% N = 864
% 0.4000    1.4600   -2.4802   -2.6944 vs -2.72\pm 0.01 from Ref. [2]
%

figure(1)
plot(md_step_ii(1:100:length(md_epot)), md_epot(1:100:length(md_epot)), 'b-') % , LineWidth=1.5
xlabel('\mbox{Time}','Interpreter','latex') % ,'fontsize',16
ylabel('$U^{\ast}$','Interpreter','latex','Rotation',1) % , 'Rotation',0
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
box on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% virial pressure (equation of state, beta*p/rho)
%
vir_pressure_long_range_corr_beta_p_over_rho = (16/3)*pi*(den/T_ref)*((2/3)*(1./rc^9) - (1./rc^3)); % from Ref. [3]. 
vir = sum(md_vir)/length(md_vir);
%
[den, T_ref, vir, 1.+ vir*(1./T_inst), 1.+ vir*(1./T_inst) + vir_pressure_long_range_corr_beta_p_over_rho]

% [den, T_ref, vir, 1.+ vir*(1/T_inst), 1.+ vir*(1/T_inst) + vir_pressure_long_range_corr_beta_p_over_rho]
%
% N = 256 
%
% 0.4000    1.4600   -0.4300    0.7057    0.4127 vs  0.41 from Ref. [2]
% 0.8800    0.9400    2.6464    3.8022    2.8010 vs  2.72 from Ref. [2]
% 0.8800    0.5910   -2.0176   -2.4202   -4.0126 vs -0.18 from Ref. [2]
% 0.8800    2.8890   12.1942    5.1040    4.7783 vs  4.36 from Ref. [2]
% 0.8500    1.2140    3.4354    3.8628    3.1141 vs  3.06 from Ref. [2]
% 0.8500    0.7600    0.7541    1.9928    0.7968 vs  0.78 from Ref. [2]
% 0.8500    0.7190    0.4564    1.6399    0.3757 vs  0.36 from Ref. [2]
% 0.8500    0.6580    0.0316    1.0474   -0.3340 vs -0.20 from Ref. [2]
% 0.8500    0.5910   -2.2757   -2.8011   -4.3392 vs -1.20 from Ref. [2]
% 0.7500    2.8490    6.8945    3.4244    3.1429 vs  3.10 from Ref. [2]
% 0.7500    1.3040    1.6204    2.2461    1.6310 vs  1.61 from Ref. [2]
% 0.7500    1.0710    0.6589    1.6143    0.8654 vs  0.89 from Ref. [2]
% 0.7500    0.8270   -0.4660    0.4339   -0.5360 vs -0.54 from Ref. [2]
% 0.6500    1.5850    1.1544    1.7428    1.3042 vs  1.25 from Ref. [2]
% 0.6500    0.9000   -0.9069   -0.0237   -0.7960 vs -0.74 from Ref. [2]
% 0.5426    1.3260   -0.2991    0.7757    0.3381 vs  0.42 from Ref. [2]
% 0.4500    1.5520   -0.2078    0.8655    0.5554 vs  0.75 from Ref. [2]
% 0.4000    1.4240   -0.4590    0.6768    0.3764 vs  0.38 from Ref. [2]
% 0.3500    1.4180   -0.5016    0.6460    0.3821 vs  0.40 from Ref. [2]


%
figure(2)
hold on
plot(md_step_ii(1:100:length(md_ave_temp)), md_ave_temp(1:100:length(md_ave_temp)), 'b', 'LineWidth',1.2)
yline(T_inst_ave, 'g--', 'LineWidth',2.5)  % computed average temperature
yline(T_ref, 'r-', 'LineWidth',2.5)        % reference temperature       
hold off
xlabel('\mbox{Time}','Interpreter','latex') % ,'fontsize',16
ylabel('\mbox{Temperature}','Interpreter','latex') % , 'Rotation',0
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
box on


%%%
return
end
%
function [x, y, z] = coord_x_y_z(A, npart, number_of_dim)
%
x = zeros(npart,1); y = zeros(npart,1); z = zeros(npart,1);
%
ijk = 0.0;
for lg = 0:1
   for i = 0:number_of_dim-1
       for j = 0:number_of_dim-1
           for k = 0:number_of_dim-1
               ijk = ijk + 1;
               x(ijk) = i * A + lg * A * 0.5;
               y(ijk) = j * A + lg * A * 0.5;
               z(ijk) = k * A;
              [ j, k, ijk, x(ijk), y(ijk), z(ijk)] ;
           end
       end
   end 
end
%
for lg = 1:2
   for i = 0:number_of_dim-1
      for j = 0:number_of_dim-1
         for k = 0:number_of_dim-1
             ijk = ijk + 1;
             x(ijk) = i * A + (2-lg) * A * 0.5;
             y(ijk) = j * A + (lg-1) * A * 0.5;
             z(ijk) = k * A + A * 0.5;
             [i, j, k, ijk, x(ijk), y(ijk), z(ijk)] ;
         end
      end
   end
end
%%%
return
end

%
function [vx, vy, vz] = init_vel(npart, T_ref, mass)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initial velocities
vx = rand(npart, 1);  %from 0 to 1
vx = 2.*vx - 1. ;    % from -1 to 1
%
vy = rand(npart, 1);
vy = 2.*vy - 1;
%
vz = rand(npart, 1);
vz = 2.*vz - 1.;
%
% enforce zero net momentum
sm_vx = sum(vx)/npart; % mean velocity 
sm_vy = sum(vy)/npart;
sm_vz = sum(vz)/npart;
%
sm_v2 = sum(vx.*vx + vy.*vy + vz.*vz);
%
vx = vx - sm_vx;
vy = vy - sm_vy;
vz = vz - sm_vz;
%
K = mass.*sm_v2/(2.*npart);   % kinetik energy per particle 
T_inst = (2/3).*K;
%
% scaling factor , fs
fs = sqrt(T_ref/T_inst);
%%%
vx = fs.*vx; 
vy = fs.*vy; 
vz = fs.*vz;

%%%
return
end

%%%
function [x, y, z] = periodic_BC(x, y, z, npart, side)
%
    % apply periodic boundary condition
    for i = 1:npart
        %
        if (x(i) < 0.);   x(i) = x(i) + side; end
        if (x(i) > side); x(i) = x(i) - side; end
        %
        if (y(i) < 0.);   y(i) = y(i) + side; end
        if (y(i) > side); y(i) = y(i) - side; end
        %
        if (z(i) < 0.);   z(i) = z(i) + side; end
        if (z(i) > side); z(i) = z(i) - side; end
        %
    end
%%%
return
end

%%%
function [ax,ay,az, epot, vir] = force_acc(x, y, z, side, sideh, rc2, npart, mass)
%
%%% internal energy, epot & accelerations, ax, ay, az
fx = zeros(npart,1); fy = zeros(npart,1); fz = zeros(npart,1);
epot = 0.;
vir = 0.;
for i = 1:npart
        %
    for j = i+1:npart
        xx = x(i) - x(j);
        yy = y(i) - y(j);
        zz = z(i) - z(j);
            %
        if (xx <-sideh); xx = xx + side; end
        if (xx > sideh); xx = xx - side; end 
            %
        if (yy <-sideh); yy = yy + side; end
        if (yy > sideh); yy = yy - side; end             
            %
        if (zz <-sideh); zz = zz + side; end
        if (zz > sideh); zz = zz - side; end 
            %
        r2 = xx * xx + yy * yy + zz * zz;
            %
        if (r2 < rc2)
            %
            r2i = 1/r2;
            r6i = r2i * r2i * r2i;
            ff  = 48 * r2i * r6i * (r6i - 0.5);
            r148 = ff;
%            virij = r2 * ff;  
            virij = 48 * (r6i * r6i - 0.5 * r6i);            
                %
             kx = xx * r148;
             fx(i) = fx(i) + kx;
             fx(j) = fx(j) - kx;
                %
             ky = yy * r148;
             fy(i) = fy(i) + ky;
             fy(j) = fy(j) - ky;    
                %
             kz = zz * r148;
             fz(i) = fz(i) + kz;
             fz(j) = fz(j) - kz;    
             %
             epot = epot + 4.*r6i * (r6i - 1.); 
             vir = vir + virij;
        end %       
    end
end
%
ax = fx./mass;
ay = fy./mass;
az = fz./mass;
%
epot = epot/npart; 
vir = vir/(3.*npart);
%%%
return
end