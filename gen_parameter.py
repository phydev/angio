# this script generate the input files 
# for a range of parameters

# which parameter would you like to vary?
alpha_p_of = True
chi_of = False

alpha_p_min, alpha_p_max, delta_alpha_p = 0.1, 3.0, 0.2
chi_min, chi_max, delta_chi = 200, 2000, 100



if(alpha_p_of):

    counters = '00'
    counteri = 0
    alpha_p = alpha_p_min
    
    while alpha_p<=alpha_p_max:

        filename = 'inp'+counters
        dir_name = 'dir'+counters

        file = open(filename, 'w')

        file.write('4.0, cell_radius \n')
        file.write('100.0, diffusion_const \n')
        file.write('1.00, interface_width \n')
        file.write('0.30, vegf_p \n')
        file.write('0.10, vegf_c \n')
        file.write('20.0, diff_oxy_length \n')
        file.write('6.25, vegf_rate \n')
        file.write('1.00, vegf_source_conc \n')
        file.write(str(alpha_p)+', prolif_rate \n')
        file.write('5.00, vessel_radius \n')
        file.write('150000, total_time_step \n')
        file.write('0.0005, dt \n')
        file.write('1000.0, chi_chemiotactic_resp \n')
        file.write('80, 80, 50, 1, 1, 1, Lx_Ly_Lz_dx_dy_dz \n')
        file.write(dir_name+', dir_name \n')
        file.write('100000, random_seed \n')
        file.write('20, number_of_boundary_points \n')
        file.write('500, source_max \n')
        file.write('0.01, vegf_grad_min \n')
        file.write('0.03, vegf_max \n')
        file.write('2.00, depletion_weight \n')
        file.write('2000, output_period \n')
        file.write('40000, extra_steps \n')
        file.write('500, max_number_of_tip_cells \n')
        file.write('T, thinning_FT \n')
        file.write('F, periodic_FT \n')
                      
        alpha_p += delta_alpha_p
        counteri += 1

        if(counteri<10):
            counters = '0'+str(counteri)
        else:
            counters = str(counteri)

        file.close()

if(chi_of):

    chi = chi_min

    while chi <= chi_max:
        chi += delta_chi
