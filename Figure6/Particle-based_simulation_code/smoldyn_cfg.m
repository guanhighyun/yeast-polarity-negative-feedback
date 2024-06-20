function smoldyn_cfg(fileprefix,directory,tstop,samplingrate,n_Cdc42,n_BemGEF,n_GAP,P3,k2b,random_seed)
% Adapted from "make_smoldyn_cfg.m" from Ramirez SA, Pablo M, Burk S,
% Lew DJ, Elston TC (2021). PLoS Comput Biol 17(7): e1008525.

% configuration file for smoldyn
cfg_name = [directory '/' fileprefix '.cfg']; 
% xyz coords of species. they'll be written to the same directory as 
% the .cfg file it is run from.
xyz_name = [fileprefix '.xyz']; 
fid=fopen(cfg_name,'w');

% Parameters for the simulation
rho = 0.05;
dt = 0.0001;
k1a = 10;
k1b = 40;
P2a = 5.3*dt;
P4a = 9.6*dt;
k4b = 40;
k5a = 36;
k5b = 13;
P7 = 256*dt;
P6a = 0.15*dt;
P6b = 110*dt;
k6c = 0.01;
Dm = 0.0025;
Dc = 15;

% Define a random number generator
fprintf(fid,'random_seed %i\n',random_seed);
% Reaction radius rho
fprintf(fid,'variable rho = %g\n',rho);
% A really small number which helps separate molecules
fprintf(fid,'variable rho_eps = 0.00001\n');
% Domain length
fprintf(fid,'variable L = 8.8623\n');

% Reaction rates
fprintf(fid,'variable k1a = %g\n',k1a);
fprintf(fid,'variable k1b = %g\n',k1b);
fprintf(fid,'variable P2a = %g\n',P2a);
fprintf(fid,'variable k2b = %g\n',k2b);
fprintf(fid,'variable P3 = %g\n',P3);
fprintf(fid,'variable P4a = %g\n',P4a);
fprintf(fid,'variable k4b = %g\n',k4b);
fprintf(fid,'variable k5a = %g\n',k5a);
fprintf(fid,'variable k5b = %g\n',k5b);
fprintf(fid,'variable P7  = %g\n',P7);
fprintf(fid,'variable P6a = %g\n',P6a);
fprintf(fid,'variable P6b = %g\n',P6b);
fprintf(fid,'variable k6c = %g\n',k6c);

% Define domain boundaries
fprintf(fid,'dim 2\n');
fprintf(fid,'boundaries x -0.2 L+0.2\n');
fprintf(fid,'boundaries y -0.2 L+0.2\n');

% Define species
% Cdc42T: Cdc42-GTP
% Cdc42Dc: Cdc42-GDP in the cytoplasm
% Cdc42Dm: Cdc42-GDP on the membrane
% BemGEFc: Bem1-GEF in the cytoplasm
% BemGEFm: Bem1-GEF on the membrane
% BemGEF42: Bem1-GEF bound with Cdc42-GTP
% GAPi: inactive GAP
% GAPa: active GAP
fprintf(fid,'species Cdc42T Cdc42Dc Cdc42Dm BemGEFc BemGEFm BemGEF42 complex_Cdc42Dm_BemGEF42 complex_Cdc42Dm_BemGEFm\n');
fprintf(fid,'species GAPa GAPi complex_GAPi_Cdc42T complex_Cdc42T_GAPa\n');

% Define diffusion rates. Dm: membrane diffusion rate. Dc: cytosolic
% diffusion rate. 
fprintf(fid,'difc Cdc42T %g\n',Dm);
fprintf(fid,'difc Cdc42Dc %g\n',Dc);
fprintf(fid,'difc Cdc42Dm %g\n',Dm);
fprintf(fid,'difc BemGEF42 %g\n',Dm);
fprintf(fid,'difc BemGEFc %g\n',Dc);
fprintf(fid,'difc BemGEFm %g\n',Dm);
fprintf(fid,'difc complex_Cdc42Dm_BemGEF42 %g\n',Dm);
fprintf(fid,'difc complex_Cdc42Dm_BemGEFm %g\n',Dm);
fprintf(fid,'difc GAPi %g\n',Dc);
fprintf(fid,'difc GAPa %g\n',Dc);
fprintf(fid,'difc complex_GAPi_Cdc42T %g\n',Dm);
fprintf(fid,'difc complex_Cdc42T_GAPa %g\n',Dm);

% Set up lists to store molecular coordinates
fprintf(fid,'molecule_lists list1 list2 list3 list4 list5 list6 list7 list8\n');
fprintf(fid,'mol_list Cdc42T list1\n');
fprintf(fid,'mol_list BemGEF42 list2\n');
fprintf(fid,'mol_list Cdc42Dm list3\n');
fprintf(fid,'mol_list Cdc42Dc list4\n');
fprintf(fid,'mol_list BemGEFm list5\n');
fprintf(fid,'mol_list BemGEFc list6\n');
fprintf(fid,'mol_list GAPa list7\n');
fprintf(fid,'mol_list GAPi list8\n');

% Define a domain. The boundary condition is periodic.
fprintf(fid,'start_surface inner_walls\n');
fprintf(fid,'action both all jump\n');
fprintf(fid,'polygon both edge\n');
fprintf(fid,'panel rect +x 0 0 L r1\n');
fprintf(fid,'panel rect -x L 0 L r2\n');
fprintf(fid,'panel rect +y 0 0 L r3\n');
fprintf(fid,'panel rect -y 0 L L r4\n');
fprintf(fid,'jump r1 front <-> r2 front\n');
fprintf(fid,'jump r3 front <-> r4 front\n');
fprintf(fid,'end_surface\n');

% Define the compartment which is required by Smoldyn to input molecules.
fprintf(fid,'start_compartment full_domain\n');
fprintf(fid,'surface inner_walls\n');
fprintf(fid,'point L/2 L/2\n');
fprintf(fid,'end_compartment\n');

% Define reactions
% Bem1-GEF associates and dissociates from the membrane
fprintf(fid,'reaction BemGEF_cmtransition BemGEFc <-> BemGEFm k1a k1b\n');

% Cdc42-GDP associates and dissociates from the membrane
fprintf(fid,'reaction Cdc42_cmtransition Cdc42Dc <-> Cdc42Dm k5a k5b\n');

% Bem1-GEFm + Cdc42Dm -> Bem1-GEFm + Cdc42T
fprintf(fid,'reaction Cdc42Dm_2_T_bindTo_BemGEFm Cdc42Dm + BemGEFm -> complex_Cdc42Dm_BemGEFm\n');
fprintf(fid,'reaction Cdc42Dm_2_T_catBy_BemGEFm complex_Cdc42Dm_BemGEFm -> Cdc42T + BemGEFm\n');
fprintf(fid,'reaction_probability Cdc42Dm_2_T_bindTo_BemGEFm P2a\n');
fprintf(fid,'binding_radius Cdc42Dm_2_T_bindTo_BemGEFm rho\n');
fprintf(fid,'reaction_probability Cdc42Dm_2_T_catBy_BemGEFm 1\n');
% Place molecules beyond their reactive radius
fprintf(fid,'product_placement Cdc42Dm_2_T_catBy_BemGEFm unbindrad rho+rho_eps\n');

% Cdc42T -> Cdc42Dm
fprintf(fid,'reaction Cdc42T_2_Cdc42Dm Cdc42T -> Cdc42Dm k2b\n');

% Cdc42T-Bem1-GEF + Cdc42Dm -> Cdc42T-Bem1-GEF + Cdc42T
fprintf(fid,'reaction Cdc42Dm_2_T_bindTo_BemGEF42 Cdc42Dm + BemGEF42 -> complex_Cdc42Dm_BemGEF42\n');
fprintf(fid,'reaction Cdc42Dm_2_T_catBy_BemGEF42 complex_Cdc42Dm_BemGEF42 -> Cdc42T + BemGEF42\n');
fprintf(fid,'reaction_probability Cdc42Dm_2_T_bindTo_BemGEF42 P3\n');
fprintf(fid,'binding_radius Cdc42Dm_2_T_bindTo_BemGEF42 rho\n');
fprintf(fid,'reaction_probability Cdc42Dm_2_T_catBy_BemGEF42 1\n');
fprintf(fid,'product_placement Cdc42Dm_2_T_catBy_BemGEF42 unbindrad rho+rho_eps\n');

% Cdc42T + Bem1-GEFm -> Cdc42T-Bem1-GEF
fprintf(fid,'reaction make_BemGEF42_fromm BemGEFm + Cdc42T <-> BemGEF42\n');
fprintf(fid,'reaction_probability make_BemGEF42_frommfwd P4a\n');
fprintf(fid,'binding_radius make_BemGEF42_frommfwd rho\n');
fprintf(fid,'reaction_rate make_BemGEF42_frommrev k4b\n');
fprintf(fid,'product_placement make_BemGEF42_frommrev unbindrad rho+rho_eps\n');

% Cdc42T + BemGEFc -> Cdc42T-Bem1-GEF
fprintf(fid,'reaction make_BemGEF42_fromc BemGEFc + Cdc42T -> BemGEF42\n');
fprintf(fid,'reaction_probability make_BemGEF42_fromc P7\n');
fprintf(fid,'binding_radius make_BemGEF42_fromc rho\n');

% Cdc42T + GAPi -> Cdc42T + GAPa
fprintf(fid,'reaction GAPi_bind_2_Cdc42T Cdc42T + GAPi -> complex_GAPi_Cdc42T\n');
fprintf(fid,'reaction_probability GAPi_bind_2_Cdc42T P6a\n');
fprintf(fid,'binding_radius GAPi_bind_2_Cdc42T rho\n');
fprintf(fid,'reaction GAPi_2_GAPa_cat_By_Cdc42T complex_GAPi_Cdc42T -> GAPa + Cdc42T\n');
fprintf(fid,'reaction_probability GAPi_2_GAPa_cat_By_Cdc42T 1\n');
fprintf(fid,'product_placement GAPi_2_GAPa_cat_By_Cdc42T unbindrad rho+rho_eps\n');

% Cdc42T + GAPa -> Cdc42Dm + GAPa
fprintf(fid,'reaction Cdc42T_bind_2_GAPa Cdc42T + GAPa -> complex_Cdc42T_GAPa\n');
fprintf(fid,'reaction_probability Cdc42T_bind_2_GAPa P6b\n');
fprintf(fid,'binding_radius Cdc42T_bind_2_GAPa rho\n');
fprintf(fid,'reaction Cdc42T_2_D_cat_By_GAPa complex_Cdc42T_GAPa -> GAPa + Cdc42Dm\n');
fprintf(fid,'reaction_probability Cdc42T_2_D_cat_By_GAPa 1\n');
fprintf(fid,'product_placement Cdc42T_2_D_cat_By_GAPa unbindrad rho+rho_eps\n');

fprintf(fid,'reaction GAPa_2_Yim GAPa -> GAPi k6c\n');

% Define the starting and ending point and time step
fprintf(fid,'time_start 0\n');
fprintf(fid,'time_stop %i\n',tstop);
fprintf(fid,'time_step %g\n',dt);

% Input initial conditions
fprintf(fid,'compartment_mol %i Cdc42Dc full_domain\n',n_Cdc42);
fprintf(fid,'compartment_mol %i BemGEFc full_domain\n',n_BemGEF);
fprintf(fid,'compartment_mol %i GAPi full_domain\n',n_GAP); 

% Suddenly reduce the positive feedback rate to a new lambda3 rate "new_lambda3"
% fprintf(fid,'cmd @ 600 set reaction_probability Cdc42Dm_2_T_bindTo_BemGEF42 %g\n',new_lambda3*dt);

% Output coordinates of molecules in the "xyz_name" every "samplingrate"
% time points.
fprintf(fid,'output_files %s\n',xyz_name);
fprintf(fid,'cmd N %i molpos Cdc42T %s\n',samplingrate,xyz_name);
fprintf(fid,'cmd N %i molpos BemGEF42 %s\n',samplingrate,xyz_name);
fprintf(fid,'cmd N %i molpos BemGEFm %s\n',samplingrate,xyz_name);
fprintf(fid,'cmd N %i molpos BemGEFc %s\n',samplingrate,xyz_name);
fprintf(fid,'cmd N %i molpos GAPa %s\n',samplingrate,xyz_name);
fprintf(fid,'cmd N %i molpos GAPi %s\n',samplingrate,xyz_name);
end