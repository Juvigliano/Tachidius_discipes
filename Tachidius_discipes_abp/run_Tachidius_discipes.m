% % close all; 
% % global pets
% % 
% % pets = {'Nitocra_spinipes_abp'};
% % check_my_pet(pets);
% % 
% % estim_options('default'); 
% % estim_options('max_step_number', 5e3); 
% % estim_options('max_fun_evals', 5e3);  
% % 
% % estim_options('pars_init_method', 2);
% % % get initial parameters from:
% % % 0 = automatized computation
% % % 1 = results_my_pet.mat file
% % % 2 = pars_init_my_pet.m file
% % estim_options('results_output', 1);
% % % 1 = show output on screen and write file
% % estim_options('method', 'no');
% % % 'no' = do not estimate
% % % 'nm' = use Nelder-Mead method
% % estim_pars;
% 
% % % write results_my_pet.mat to pars_init_my_pet.m
% % mat2pars_init('Nitocra_spinipes_abp')   
close all;
global pets

pets = {'Tachidius_discipes'};
check_my_pet(pets);

estim_options('default');
estim_options('max_step_number',2e3);
estim_options('max_fun_evals',5e3);

estim_options('pars_init_method', 2);
estim_options('results_output', 5);
estim_options('method', 'no');

estim_pars;
% mat2pars_init;

