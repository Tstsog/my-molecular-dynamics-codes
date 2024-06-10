% This Matlab code performs the molecular dynamics simulation. The code is
% my own written Matlab version of the Microcanonical Molecular-Dynamics
% Program written in Fortran in Ref. [1]. Results from this calculation
% are compared with those obtained in Ref. [1]. 
%
% An original code in Fortran is Program Listing PL 1. Microcanonical Molecular-Dynamics Program in 
% the Appendix of Ref.[1]. 
%
% Ref. [1] D. Heermann, "Computer Simulation Methods in Theoretical
% Physics", 2nd edition, (1989)
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
%
% June 5, 2024 & University of North Dakota 
%
%%%%
function [] = heermann_md_code
clc; clear; 
format long
%
dim_parm = 4.;             % dimensional parameter
npart = 4 * dim_parm^3;    % number of particle 
den = 0.83134;             % number density 
side = (npart/den)^(1/3);  % box's size length
tref = 0.722;              % reduced reference temperature  
rcoff = 2.5;               % cut-off parameter in the Lennard-Jones potential 
h = 0.064;                 % basis time step 
%
A = side/dim_parm;
sideh = side * 0.5;
hsq = h * h;
hsq2 = hsq * 0.5;
rcoffs = rcoff * rcoff;
tscale = 16.0 / (1.0 * npart - 1.0);
vaver = 1.13 * sqrt(tref / 24.0);
%
x = zeros(npart, 1); y = zeros(npart, 1); z = zeros(npart, 1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
% SET UP FCC LATTICE FOR THE ATOMS INSIDE THE BOX
ijk = 0.;
for lg = 0:1
    for i = 0:dim_parm-1
        for j = 0:dim_parm-1
            for k = 0:dim_parm-1
                ijk = ijk + 1;
                x(ijk) = i * A + lg * A * 0.5;
                y(ijk) = j * A + lg * A * 0.5;
                z(ijk) = k * A;
%                [i,j,k, y(ijk), z(ijk)];
            end
        end
    end
end
%
for lg = 1:2
    for i = 0:dim_parm-1
        for j = 0:dim_parm-1
            for k = 0:dim_parm-1
                ijk = ijk + 1;
                x(ijk) = i * A + (2-lg) * A * 0.5;
                y(ijk) = j * A + (lg-1) * A * 0.5;
                z(ijk) = k * A + A * 0.5;
%                [i,j,k, y(ijk), z(ijk)];
            end
        end
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%   ASSIGN VELOCITIES DISTRIBUTED NORMALLY
[vh_x, vh_y, vh_z] = mxwell(npart, tscale, tref, h);

%%% 
f_x = zeros(npart,1); f_y = zeros(npart,1); f_z = zeros(npart,1);

fileID_save_data_1 = fopen('heermann_md_code.txt','w');

md_step = 1000;  % number of MD step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF THE ACTUAL MOLECULAR DYNAMICS PROGRAM.
% THE EQUATIONS OF MOTION ARE INTEGRATED USING THE 'SUMMED FORM' 
% (G. DAHLQUIST AND A. BJOERK, NUMERICAL METHODS, PRENTICE HALL 
% ( 1 9 7 4 » . EVERY lREP'TH STEP THE VELOCITIES ARE RESCALED SO AS
% TO GIVE THE SPECIFIED TEMPERATURE !REF.
% VERSION 1 . 0 AUGUST 1986 DIETER W. HEERMANN

for ii = 1:md_step

%%% ADVANCE POSITIONS ONE BASIC TIME STEP
for i = 1:npart 
    x(i) = x(i) + vh_x(i) + f_x(i);
    y(i) = y(i) + vh_y(i) + f_y(i);    
    z(i) = z(i) + vh_z(i) + f_z(i);        
    [i, x(i), y(i), z(i)];
end
%
%%% APPLY PERIODIC BOUNDARY CONDITIONS
for i = 1:npart 
    %
    if ( x(i) < 0.); x(i) = x(i) + side; end
    if ( x(i) > side); x(i) = x(i) - side; end
    %
    if ( y(i) < 0.); y(i) = y(i) + side; end
    if ( y(i) > side); y(i) = y(i) - side; end
    %
    if ( z(i) < 0.); z(i) = z(i) + side; end
    if ( z(i) > side); z(i) = z(i) - side; end
    %
    [i, x(i), y(i), z(i)];
end
%
%%% OMPUTE THE PARTIAL VELOCITIES
for i = 1:npart
    vh_x(i) = vh_x(i) + f_x(i);
    vh_y(i) = vh_y(i) + f_y(i);
    vh_z(i) = vh_z(i) + f_z(i);  
    %
    [i, vh_x(i), vh_y(i), vh_z(i)];
end

vir = 0.;
epot = 0.;
%
for i = 1:npart
    for j = i+1:npart
        xx = x(i) - x(j);
        yy = y(i) - y(j);        
        zz = z(i) - z(j);                
        %
        %%% MINIMUM IMAGE CONVENTION 
        if (xx < -sideh); xx = xx + side; end 
        if (xx >  sideh); xx = xx - side; end     
        %
        if (yy < -sideh); yy = yy + side; end 
        if (yy >  sideh); yy = yy - side; end             
        %
        if (zz < -sideh); zz = zz + side; end 
        if (zz >  sideh); zz = zz - side; end     
        %
        rd = xx * xx + yy * yy + zz * zz;
        %
        if (rd < rcoffs)
            epot = epot + rd^(-6) - rd^(-3);
            r148 = rd^(-7) - 0.5 * rd^(-4);
            vir = vir - rd * r148;
            %
            kx = xx * r148;
            f_x(i) = f_x(i) + kx;
            f_x(j) = f_x(j) - kx;      
            %
            ky = yy * r148;
            f_y(i) = f_y(i) + ky;
            f_y(j) = f_y(j) - ky;      
            %
            kz = zz * r148;
            f_z(i) = f_z(i) + kz;
            f_z(j) = f_z(j) - kz;                  

        end
        %
    end
end
%
for i = 1:npart
    f_x(i) = f_x(i) * hsq2;
    f_y(i) = f_y(i) * hsq2;    
    f_z(i) = f_z(i) * hsq2;    
end
%
%%% COMPUTE THE VELOCITIES 
for i = 1:npart
    vh_x(i) = vh_x(i) + f_x(i);
    vh_y(i) = vh_y(i) + f_y(i);
    vh_z(i) = vh_z(i) + f_z(i);  
    %
    [i, vh_x(i), vh_y(i), vh_z(i)];
end
%
%%% COMPUTE THE KINETIC ENERGY 
ekin = 0.;
for i = 1:npart
    ekin = ekin + vh_x(i)*vh_x(i) + vh_y(i)*vh_y(i) + vh_z(i)*vh_z(i);
end
ekin = ekin/hsq ;
%
%%% COMPUTE THE AVERAGE VELOCITY 
vel = 0.;
count  = 0.;
sq = zeros(npart,1); sqt = zeros(npart,1);
for i = 1:npart
    sq(i) = sqrt(vh_x(i)*vh_x(i) + vh_y(i)*vh_y(i) + vh_z(i)*vh_z(i));
    sqt(i) = sq(i)/h;
    if (sqt(i) > vaver)
        count = count + 1;
    end
    vel = vel + sq(i);
end
vel = vel/h;
%
% scaling 
ts = tscale * ekin;
sc = tref/ts;
sc = sqrt(sc);
%
for i = 1:npart 
    vh_x(i) = vh_x(i) * sc;
    vh_y(i) = vh_y(i) * sc;    
    vh_z(i) = vh_z(i) * sc;    
%    [i, vh_x(i), vh_y(i), vh_z(i)];
end
ekin = tref/tscale; 

%
[ekin, epot, vir];
ek = 24 * ekin;
epot = 4 * epot;
etot = ek + epot;
temp = tscale * ekin;
press = den*16*(ekin - vir)/npart;
vel = vel/npart; 
rp = (count/npart) * 100;

output = [ii, ek, epot, etot, temp, press, vel, rp ];

fprintf(fileID_save_data_1, '%4.6f \t %4.8f \t %4.8f \t %4.8f \t %4.8f \t %4.6f \t %4.6f \t %8.6f\n', output); 

end 
fclose(fileID_save_data_1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_md = fopen('heermann_md_code.txt','r');
data_read_md = textscan(data_md, '%f %f %f %f %f %f %f %f') ;
md_step_ii = data_read_md{1};
md_pot_kin = data_read_md{2};
md_pot_en = data_read_md{3};
md_tot_en = data_read_md{4};
md_vir_press = data_read_md{6};
md_ave_vel = data_read_md{7};
md_rp = data_read_md{8};

%%%
figure(1)
plot(md_step_ii, md_pot_en, 'b-', LineWidth=1.5)
xlabel('\mbox{MD steps}','Interpreter','latex') % ,'fontsize',16
ylabel('$U^{\ast}$','Interpreter','latex', 'Rotation',1) %
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',20)
box on
%
ave_kin_en = sum(md_pot_kin)/length(md_pot_kin) 
ave_pot_en = sum(md_pot_en)/length(md_pot_en);
ave_tot_en = sum(md_tot_en)/length(md_tot_en);
ave_press = sum(md_vir_press)/length(md_ave_vel);
ave_vel = sum(md_ave_vel)/length(md_ave_vel);
ave_rp = sum(md_rp)/length(md_rp) 
%
%
[ave_kin_en] % 2.761650000000004e+02    vs   279.13 from Ref.[1]  
[ave_pot_en] % -1.431099468893100e+03   vs -1421.98 from Ref.[1]
[ave_tot_en] % -1.154934468893101e+03   vs -1142.92 from Ref.[1]
[ave_press]  % 0.518615560000000        vs 
[ave_vel]    % 0.195136427000000        vs 0.1965 from Ref.[1]
[ave_rp]     % 46.376171874999997       vs 47.08 from Ref.[1]

%%%
return
end

%%%

function [vh_x, vh_y, vh_z] = mxwell(npart, tscale, tref, h)
%
u1_rnd_x = rand(npart,1);
u1_rnd_y = rand(npart,1);
u1_rnd_z = rand(npart,1);
%
vh_x = 2.*u1_rnd_x - 1.;
vh_y = 2.*u1_rnd_y - 1.;
vh_z = 2.*u1_rnd_z - 1.;

%
ekin_x = 0.;
sp_x = 0.;
for i = 1:npart
    sp_x = sp_x + vh_x(i);
end
sp_x = sp_x/npart;
for i = 1:npart
    vh_x(i) = vh_x(i) - sp_x;
    ekin_x = ekin_x + vh_x(i) * vh_x(i);
end
%ekin_x;
%
sp_y = 0.;
for i = 1:npart
    sp_y = sp_y + vh_y(i);
end
sp_y = sp_y/npart;
ekin_y = 0.;
for i = 1:npart
    vh_y(i) = vh_y(i) - sp_y;
    ekin_y = ekin_y + vh_y(i) * vh_y(i);
end
%ekin_y;
%
sp_z = 0.;
for i = 1:npart
    sp_z = sp_z + vh_z(i);
end
sp_z = sp_z/npart;
ekin_z = 0.;
for i = 1:npart
    vh_z(i) = vh_z(i) - sp_z;
    ekin_z = ekin_z + vh_z(i) * vh_z(i);
end
%
ekin = ekin_x + ekin_y + ekin_z;
%
ts = tscale * ekin;
sc = tref/ts;
sc = sqrt(sc);
sc = sc * h;
for i = 1:npart 
    vh_x(i) = vh_x(i) * sc;
    vh_y(i) = vh_y(i) * sc;    
    vh_z(i) = vh_z(i) * sc;    
%    [i, vh_x(i), vh_y(i), vh_z(i)];
end
%%%

return
end
