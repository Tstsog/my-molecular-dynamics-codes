% This Matlab code performs the molecular dynamics (MD) simulation for one-component Lennard-Jones (LJ) fluid. 
% A velocity-Verlet scheme is used in numerical integration. 
% Results for a potential (pot) energy (en) per particle from this calculation
% are compared with those obtained by Verlet, Ref. [2]. 
%
% This matlab code is written in similar way that is used in Fortran code in 
% the Appendix of Ref.[1]. 
%
% Ref. [1] D. Heermann, "Computer Simulation Methods in Theoretical Physics", 2nd edition, (1989);
% Ref. [2] L. Verlet, Phys. Rev. v159, p98 (1967); 
% Ref. [3] D. Frenkel and B. Smit, "Understanding Moleculer Simulation", Acedmic Press (2002);
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
%
% July 2, 2024 & University of North Dakota 
%
function [] = classical_MD_for_LJ_fluid
clc; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
number_of_dim = 4;           % numerical parameter & you may change it
npart = 4 * number_of_dim^3; % number of particle 
den = 0.45;    %  reduced density 
T_ref = 1.744; %  reduced reference temperature
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
    [ax1, ay1, az1, epot] = force_acc(x, y, z, side, sideh, rc2, npart, mass); % velocity Verlet scheme: stage 2
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
    output = [ii*dt, epot];

    fprintf(fileID_save_data_1, '%4.4f \t %8.4f\n', output); 

end
fclose(fileID_save_data_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_md_data = fopen('classical_MD_for_LJ_fluid.txt', 'r');               % 
read_md_data = textscan(read_md_data, '%f %f');
md_step_ii = read_md_data{1};
md_epot = read_md_data{2};

%%%
Potential_energy_tail_correction_per_particle = ((8/3)*pi*den)*((1/3)*(1./rc^9) - (1./rc^3));   % from Ref. [3]. 
Potential_energy_per_partcile_without_tail_correction = sum(md_epot)/length(md_epot);
Potential_energy_per_partcile_with_tail_correction = sum(md_epot)/length(md_epot) + Potential_energy_tail_correction_per_particle;
%%%
[den, T_ref, Potential_energy_per_partcile_without_tail_correction, Potential_energy_per_partcile_with_tail_correction]

%[den, T_ref, Potential_energy_per_partcile_without_tail_correction, Potential_energy_per_partcile_with_tail_correction]
%
% 0.5426    1.4040   -3.3329   -3.6234 vs -3.63 from Ref. [2]
% 0.8500    2.2020   -4.2777   -4.7329 vs -4.76 from Ref. [2]
% 0.8800    0.9400   -5.5625   -6.0336 vs -5.84 from Ref. [2]
% 0.4500    1.7440   -2.6625   -2.9035 vs -2.90 from Ref. [2]

figure(1)
plot(md_step_ii, md_epot, 'b-') % , LineWidth=1.5
xlabel('\mbox{Time}','Interpreter','latex') % ,'fontsize',16
ylabel('$U^{\ast}$','Interpreter','latex','Rotation',1) % , 'Rotation',0
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',18)
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
function [ax,ay,az, epot] = force_acc(x, y, z, side, sideh, rc2, npart, mass)
%
%%% internal (potential) energy, epot & accelerations, ax, ay, az
fx = zeros(npart,1); fy = zeros(npart,1); fz = zeros(npart,1);
epot = 0.;
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
        end %       
    end
end
%
ax = fx./mass;
ay = fy./mass;
az = fz./mass;
%
epot = epot/npart; 
%%%
return
end