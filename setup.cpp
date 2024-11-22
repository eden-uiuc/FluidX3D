
	#include "setup.hpp"
	#include "info.hpp"

	void main_setup() { // Airfoil Test; required extensions in defines.hpp: D2Q9, FP16C, EQUILIBRIUM_BOUNDARIES, INTERACTIVE_GRAPHICS, SUBGRID, FORCE_FIELD
    // ################################################################## define simulation box size, viscosity and volume force ###################################################################

    // LBM Parameters
    const uint memory = 20000u;
    const uint T_Steps = 15u;
    const uint x_length = 2048u;
    const uint y_length = 1024u;
    const uint z_length = 1u;

    //Flow Parameters
    const float Re = 6e6f;                                    //Reynolds Number
    const float si_u = 51.45f;                               //Velocity, Initial Input was 88.8f
    const float si_rho = 1.225f;                            //Density of Fluid

    //Airfoil Parameters
    const float c = float(x_length)/20.0f;               //equivalent to si_x
    const float t = 0.1f * c;                   //thickness
    const float p = 0.6f;                //position of maximum camber
    const float m = 0.05f;                        //maximum camber
    const float aoa = 16.3759f;                      //angle of attack

    //convergence checks
    float Cl_pen = 0.0f;
    float Cl_last = 0.0f;
    float cd_pen = 0.0f;
    float cd_last = 0.0f;

    //Setup Units – Convert SI to LBM Units
    units.set_m_kg_s(c, 0.05f, 1.0f, 1.0f, si_u, si_rho);


    LBM lbm(x_length, y_length, z_length, units.nu_from_Re(Re, c, si_u));
    // ###################################################################################### define geometry ######################################################################################
    const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
        if(airfoil_cam(x, y, z, float3(Nx/2u, Ny/2u, Nz/2u), c, t, p, m, aoa)) {lbm.flags[n] = TYPE_S;}
        else lbm.u.x[n] = 0.05f;
        if(x==0u||x==Nx-1u||y==0u||y==Ny-1u) lbm.flags[n] = TYPE_E; // all non periodic
    }); // ####################################################################### run simulation, export images and data ##########################################################################
    lbm.graphics.visualization_modes = VIS_FIELD;
    lbm.graphics.slice_mode = 3;
    lbm.run(units.t(0.1));
    const string path = get_exe_path()+"FP16C/"+to_string(memory)+"MB/";
    lbm.write_status(path);
    write_file(path+"Cd"+to_string(aoa)+".dat", "# t	Cd");
    write_file(path+"Cl"+to_string(aoa)+".dat", "# t Cl");

    while(lbm.get_t()<=units.t(T_Steps)) { // main simulation loop

        Clock clock;
        lbm.calculate_force_on_boundaries();
        lbm.F.read_from_device();

        const float3 lbm_force = lbm.calculate_force_on_object(TYPE_S);
        const float force_y = units.si_F(lbm_force.y);
        const float force_x = units.si_F(lbm_force.x);

        const float Cl = force_x / (0.5f * si_rho * sq(si_u));          // calculation of lift coefficient; Multiplication by Cross-Section Area Missing?
        const float Cd = force_y / (0.5f * si_rho * sq(si_u));          // expect Cd to be too large by a factor 1.3-2.0x; need wall model

        write_line(path + "Cd"+to_string(aoa)+".dat", to_string(lbm.get_t()) + "	" + to_string(Cd, 6u) + "");
        write_line(path + "Cl"+to_string(aoa)+".dat", to_string(lbm.get_t()) + "	" + to_string(Cl, 6u) + "");

        /*
        if(lbm.graphics.next_frame(T_Steps, 25.0f)) { // render enough frames for 25 seconds of 60fps video
            lbm.graphics.set_camera_centered(-90.0f, 90.0f, 60.0f, 1.750f); // set camera position
            lbm.graphics.write_frame(get_exe_path()+"export/camera_1/"); // export image from camera position
        }
        */

        const float threshold = 1E-4f;

        lbm.run(1u, T_Steps); // run 1 LBM time step
        if (converged(Cl_pen, Cl_last, Cl, threshold) && converged(cd_pen, cd_last, Cd, threshold))
        {
            exit(0);
        }
        else
        {
            Cl_last = Cl;
            Cl_pen = Cl_last;
            cd_last = Cd;
            cd_pen = cd_last;
        }
    }
}

