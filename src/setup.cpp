#include "setup.hpp"
#include "info.hpp"

void main_setup() { // Airfoil Test; required extensions in defines.hpp: D2Q9, FP16C, EQUILIBRIUM_BOUNDARIES, INTERACTIVE_GRAPHICS, SUBGRID, FORCE_FIELD
    // ################################################################## define simulation box size, viscosity and volume force ###################################################################

    // Geometry Parameters
    const uint XY_aspect_ratio = 2u;
    const uint x_length = 2048u;
    const uint y_length = x_length/XY_aspect_ratio;
    const uint z_length = 1u;

    // Flow Parameters
    const float si_Re = 6e6f;                                             //Reynolds Number
    const float si_mu = 1.852e-5f;                                        //Viscosity of Air at 20 C
    const float si_rho = 1.225f;                 						  //Density of Fluid
    const float si_u = units.si_u_from_si_Re(si_Re, 1.0, si_mu, si_rho);  //Calculate velocity from Reynold's Number

    // Airfoil Parameters
    const float c = float(x_length)/20.0f;       //equivalent to si_x
    const float t = 0.12f * c;                   //thickness
    const float p = 0.0f;                		 //position of maximum camber
    const float m = 0.0f;                        //maximum camber
    const float aoa = 7.96346f;                  //angle of attack

    // Setup Units – Convert SI to LBM Units
    units.set_m_kg_s(c, 0.05f, 1.0f, 1.0f, si_u, si_rho);
    const ulong lbm_T = units.t(5.0f);

	// Initialize LBM using units construct to calculate LBM kinematic viscosity
    LBM lbm(x_length, y_length, z_length, units.nu(si_mu/si_rho));
    // ###################################################################################### define geometry ######################################################################################
    const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
        if(airfoil_cam(x, y, z, float3(Nx/2u, Ny/2u, Nz/2u), c, t, p, m, aoa)) {lbm.flags[n] = TYPE_S;}
        else lbm.u.x[n] = 0.05f;
        if(x==0u||x==Nx-1u||y==0u||y==Ny-1u) lbm.flags[n] = TYPE_E; // all non periodic
    }); // ####################################################################### run simulation, export images and data ##########################################################################

//    lbm.graphics.visualization_modes = VIS_FIELD;
//    lbm.graphics.slice_mode = 3;

    lbm.run(0u);
    const string path = get_exe_path()+"FP16C/NACA_0012/";
    lbm.write_status(path);
    write_file(path+"Cd.dat", "# t\tCd\n");
    write_file(path+"Cl.dat", "# t\tCl\n");

    float Cl_smooth_old = 0.0f;
    float Cd_smooth_old = 0.0f;
    const float smoothing = 5.0f;

    bool done = false;
    const float tol = 1e-8f;

    while(lbm.get_t()<=lbm_T) { // main simulation loop
        Clock clock;
        lbm.calculate_force_on_boundaries();
        lbm.F.read_from_device();
        const float3 lbm_force = lbm.calculate_force_on_object(TYPE_S);
        const float force_y = units.si_F(lbm_force.y);
        const float force_x = units.si_F(lbm_force.x);
        const float q = 0.5f * si_rho * sq(si_u);
        const float Cl = force_y / (q * units.si_x(c));          // calculation of lift coefficient; Multiplication by Cross-Section Area Missing?
        const float Cd = force_x / (q * units.si_x(c));          // expect Cd to be too large by a factor 1.3-2.0x; need wall model

        float Cl_smooth = Cl * (smoothing / (1 + lbm.get_t())) + Cl_smooth_old * (1 - smoothing / (1 + lbm.get_t()));
        float Cd_smooth = Cd * (smoothing / (1 + lbm.get_t())) + Cd_smooth_old * (1 - smoothing / (1 + lbm.get_t()));

        if (fabs(Cl_smooth - Cl_smooth_old) < tol && fabs(Cd_smooth - Cd_smooth_old) < tol){
            write_line(path + "Cd.dat", to_string(lbm.get_t()) + "\t" + to_string(Cd, 6u) + "\n");
            write_line(path + "Cl.dat", to_string(lbm.get_t()) + "\t" + to_string(Cl, 6u) + "\n");
            done = true;
        }
        if (done){
            exit(0);
        } else {
            Cl_smooth_old = Cl_smooth;
            Cd_smooth_old = Cd_smooth;
        }

//        if(lbm.graphics.next_frame(lbm_T, 5.0f)) { // render enough frames for 25 seconds of 60fps video
//            lbm.graphics.set_camera_centered(-90.0f, 90.0f, 60.0f, 1.750f); // set camera position
//            lbm.graphics.write_frame(get_exe_path()+"export/camera_1/"); // export image from camera position
//        }
        lbm.run(1u, lbm_T); // run 1 LBM time step
    }
}

