#include "setup.hpp"
#include "info.hpp"

void main_setup() { // Airfoil Test; required extensions in defines.hpp: D2Q9, FP16C, EQUILIBRIUM_BOUNDARIES, INTERACTIVE_GRAPHICS, SUBGRID, FORCE_FIELD
    // ################################################################## define simulation box size, viscosity and volume force ###################################################################
    const uint memory = 200u;
    const uint T_Steps = 15000u;
    const uint x_length = 1024u;
    const uint y_length = 2048u;
    const uint z_length = 1u;
    const float c = float(y_length)/5.0f;
    const float t = 0.10f * c;
    const float Re = 250.0f;
    const float si_u = 0.10f;
    const float si_rho = 1.0f;
    const float p = 0.40f;
    const float m = 0.06f;
    //const float aoa = 45.0f;
    const float si_A = x_length*y_length;
    LBM lbm(x_length, y_length, z_length, units.nu_from_Re(Re, x_length, si_u));
    // ###################################################################################### define geometry ######################################################################################
    const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
        if(airfoil_cam(x, y, z, float3(Nx/2u, Ny/2u, Nz/2u), c, t, p, m)) {lbm.flags[n] = TYPE_S;}
        else lbm.u.y[n] = si_u;
        if(x==0u||x==Nx-1u||y==0u||y==Ny-1u) lbm.flags[n] = TYPE_E; // all non periodic
    }); // ####################################################################### run simulation, export images and data ##########################################################################
    lbm.graphics.visualization_modes = VIS_FIELD;
    lbm.graphics.slice_mode = 3;
    lbm.run();
    const string path = get_exe_path()+"FP16C/"+to_string(memory)+"MB/";
    lbm.write_status(path);
    write_file(path+"Cd.dat", "# t\tCd\n");
    while(lbm.get_t()<=units.t(T_Steps)) { // main simulation loop
        Clock clock;
        lbm.calculate_force_on_boundaries();
        const float3 lbm_force = lbm.calculate_force_on_object(TYPE_S);
        const float force = units.si_F(lbm_force.x);
        const float Cd = units.si_F(lbm_force.x) / (0.5f * si_rho * sq(si_u) * si_A); // expect Cd to be too large by a factor 1.3-2.0x; need wall model                                                                         ");
        write_line(path + "Cd.dat", to_string(lbm.get_t()) + "\t" + to_string(Cd, 3u) + "\n");
    }
    lbm.flags.write_device_to_vtk();
    lbm.write_status(path);
}