void export_containers_methods(jlcxx::Module& mod)
{
    using namespace voro;

    // container_base methods
    
    auto __cxxwrap_bounds = [] (auto &con)
    {
        return std::make_tuple(con.ax, con.ay, con.az, con.bx, con.by, con.bz);
    };

    auto __cxxwrap_bounds_triclinic = [] (auto &con)
    {
        return std::make_tuple(con.bx, con.bxy, con.by, con.bxz, con.byz, con.bz);
    };

    auto __cxxwrap_periodic = [] (auto &con)
    {
        return std::make_tuple(con.xperiodic, con.yperiodic, con.zperiodic);
    };

    auto __cxxwrap_point_inside = [] (auto &w, double x, double y, double z) 
    {
        return w.point_inside(x, y, z);
    };

    auto __cxxwrap_draw_domain_gnuplot_file = [] (auto& con, FILE* fp) { con.draw_domain_gnuplot(fp); };

    auto __cxxwrap_draw_domain_gnuplot_filename = [] (auto& con, const char* filename) { con.draw_domain_gnuplot(filename); };

    auto __cxxwrap_draw_domain_pov_file = [] (auto& con, FILE* fp) { con.draw_domain_pov(fp); };

    auto __cxxwrap_draw_domain_pov_filename = [] (auto& con, const char* filename) { con.draw_domain_pov(filename); };

    auto __cxxwrap_total_particles = [] (auto &con) 
    {
        return con.total_particles();
    };

    auto __cxxwrap_clear = [] (auto &con) 
    {
        return con.clear();
    };

    auto __cxxwrap_import_file = [] (auto& con, FILE* fp)
    {
        con.import(fp);
    };

    auto __cxxwrap_ordered_import_file = [] (auto& con, particle_order& ord, FILE* fp)
    {
        con.import(ord, fp);
    };

    auto __cxxwrap_import_filename = [] (auto& con, const char* filename)
    {
        con.import(filename);
    };

    auto __cxxwrap_ordered_import_filename = [] (auto& con, particle_order& ord, const char* filename)
    {
        con.import(ord, filename);
    };
    
    auto __cxxwrap_compute_all_cells = [] (auto& con) { con.compute_all_cells(); };

    auto __cxxwrap_sum_cell_volumes = [] (auto& con) { con.sum_cell_volumes(); };

    auto __cxxwrap_draw_particles_file = [] (auto& con, FILE* fp) { con.draw_particles(fp); };

    auto __cxxwrap_draw_particles_filename = [] (auto& con, const char* filename) { con.draw_particles(filename); };

    auto __cxxwrap_draw_particles_pov_file = [] (auto& con, FILE* fp) { con.draw_particles_pov(fp); };

    auto __cxxwrap_draw_particles_pov_filename = [] (auto& con, const char* filename) { con.draw_particles_pov(filename); };

    auto __cxxwrap_draw_cells_file = [] (auto& con, FILE* fp) { con.draw_cells_gnuplot(fp); };

    auto __cxxwrap_draw_cells_filename = [] (auto& con, const char* filename) { con.draw_cells_gnuplot(filename); };

    auto __cxxwrap_draw_cells_pov_file = [] (auto& con, FILE* fp) { con.draw_cells_pov(fp); };

    auto __cxxwrap_draw_cells_pov_filename = [] (auto& con, const char* filename) { con.draw_cells_pov(filename); };

    auto __cxxwrap_compute_cell = [] (voronoicell_neighbor& vc, auto& con, auto& vl)
    {
        return con.compute_cell(vc, vl);
    };

    auto __cxxwrap_find_voronoi_cell = [](auto& con, double x, double y, double z)
    {
        int pid;
        double rx, ry, rz;
        bool found = con.find_voronoi_cell(x, y, z, rx, ry, rz, pid);
        return std::make_tuple(rx, ry, rz, pid, found);
    };

    auto __cxxwrap_compute_ghost_cell = [](voronoicell_neighbor &vc, auto& con, double x, double y, double z)
    {
        return con.compute_ghost_cell(vc, x, y, z);
    };
    
    auto __cxxwrap_compute_ghost_cell_poly = [](voronoicell_neighbor &vc, auto& con, double x, double y, double z, double r)
    {
        return con.compute_ghost_cell(vc, x, y, z, r);
    };

    auto __cxxwrap_add_wall = [] (auto& con, auto& wall) { con.add_wall(wall); };

    auto __cxxwrap_point_inside_walls = [] (auto &con, double x, double y, double z) 
    {
        return con.point_inside_walls(x, y, z);
    };

    auto __cxxwrap_apply_walls = [] (voronoicell_neighbor& vc, auto& con, double x, double y, double z) {
        return con.apply_walls(vc, x, y, z);
    };

    // Class Container
    mod.method(
        "__cxxwrap_bounds",
        static_cast<std::tuple<double,double,double,double,double,double> (*)(container&)>(
            __cxxwrap_bounds
        )
    );
    mod.method(
        "__cxxwrap_periodic",
        static_cast<std::tuple<bool,bool,bool> (*)(container&)>(
            __cxxwrap_periodic
        )
    );
    mod.method(
        "__cxxwrap_point_inside",
        static_cast<bool (*)(container&, double, double, double)>(
            __cxxwrap_point_inside
        )
    );
    mod.method(
        "__cxxwrap_draw_domain_gnuplot",
        static_cast<void (*)(container&, FILE*)>(
            __cxxwrap_draw_domain_gnuplot_file
        )
    );
    mod.method(
        "__cxxwrap_draw_domain_gnuplot",
        static_cast<void (*)(container&, const char*)>(
            __cxxwrap_draw_domain_gnuplot_filename
        )
    );
    mod.method(
        "__cxxwrap_draw_domain_pov",
        static_cast<void (*)(container&, FILE*)>(
            __cxxwrap_draw_domain_pov_file
        )
    );
    mod.method(
        "__cxxwrap_draw_domain_pov",
        static_cast<void (*)(container&, const char*)>(
            __cxxwrap_draw_domain_pov_filename
        )
    );
    mod.method(
        "__cxxwrap_total_particles",
        static_cast<int (*)(container&)>(
            __cxxwrap_total_particles
        )
    );
    mod.method(
        "__cxxwrap_clear!",
        static_cast<void (*)(container&)>(
            __cxxwrap_clear
        )
    );
    mod.method(
        "__cxxwrap_put!",
        static_cast<void (*)(container&, int, double, double, double)>(
            [](container& con, int i, double x, double y, double z)
            {
                con.put(i, x, y, z);
            }
        )
    );
    mod.method(
        "__cxxwrap_put!",
        static_cast<void (*)(container&, particle_order&, int, double, double, double)>(
            [](container& con, particle_order& ord, int i, double x, double y, double z)
            {
                con.put(ord, i, x, y, z);
            }
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container&, FILE*)>(
            __cxxwrap_import_file
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container&, particle_order&, FILE*)>(
            __cxxwrap_ordered_import_file
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container&, const char*)>(
            __cxxwrap_import_filename
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container&, particle_order&, const char*)>(
            __cxxwrap_ordered_import_filename
        )
    );
    mod.method(
        "__cxxwrap_compute_all_cells",
        static_cast<void (*)(container&)>(
            __cxxwrap_compute_all_cells
        )
    );
    mod.method(
        "__cxxwrap_sum_cell_volumes",
        static_cast<void (*)(container&)>(
            __cxxwrap_sum_cell_volumes
        )
    );
    mod.method(
        "__cxxwrap_draw_particles",
        static_cast<void (*)(container&, FILE*)>(
            __cxxwrap_draw_particles_file
        )
    );
    mod.method(
        "__cxxwrap_draw_particles",
        static_cast<void (*)(container&, const char*)>(
            __cxxwrap_draw_particles_filename
        )
    );
    mod.method(
        "__cxxwrap_draw_particles_pov",
        static_cast<void (*)(container&, FILE*)>(
            __cxxwrap_draw_particles_pov_file
        )
    );
    mod.method(
        "__cxxwrap_draw_particles_pov",
        static_cast<void (*)(container&, const char*)>(
            __cxxwrap_draw_particles_pov_filename
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_gnuplot",
        static_cast<void (*)(container&, FILE*)>(
            __cxxwrap_draw_cells_file
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_gnuplot",
        static_cast<void (*)(container&, const char*)>(
            __cxxwrap_draw_cells_filename
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_pov",
        static_cast<void (*)(container&, FILE*)>(
            __cxxwrap_draw_cells_pov_file
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_pov",
        static_cast<void (*)(container&, const char*)>(
            __cxxwrap_draw_cells_pov_filename
        )
    );
    mod.method(
        "__cxxwrap_find_voronoi_cell",
        static_cast<std::tuple<double, double, double, int, bool> (*)(container&, double, double, double)>(
            __cxxwrap_find_voronoi_cell
        )
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container&, c_loop_all&)>(
            __cxxwrap_compute_cell
        )
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container&, c_loop_order&)>(
            __cxxwrap_compute_cell
        )
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container&, c_loop_subset&)>(
            __cxxwrap_compute_cell
        )
    );
    mod.method(
        "__cxxwrap_point_inside_walls",
        static_cast<bool (*)(container&, double, double, double)>(
            __cxxwrap_point_inside_walls
        )
    );
    mod.method(
        "__cxxwrap_add_wall!",
        static_cast<void (*)(container&, wall_plane&)>(
            __cxxwrap_add_wall
        )
    );
    mod.method(
        "__cxxwrap_add_wall!",
        static_cast<void (*)(container&, wall_sphere&)>(
            __cxxwrap_add_wall
        )
    );
    mod.method(
        "__cxxwrap_add_wall!",
        static_cast<void (*)(container&, wall_cylinder&)>(
            __cxxwrap_add_wall
        )
    );
    mod.method(
        "__cxxwrap_add_wall!",
        static_cast<void (*)(container&, wall_cone&)>(
            __cxxwrap_add_wall
        )
    );
    mod.method(
        "__cxxwrap_apply_walls!",
        static_cast<bool (*)(voronoicell_neighbor&, container&, double, double, double)>(
            __cxxwrap_apply_walls
        )
    );

    // Class Container Poly
    mod.method(
        "__cxxwrap_bounds",
        static_cast<std::tuple<double,double,double,double,double,double> (*)(container_poly&)>(
            __cxxwrap_bounds
        )
    );
    mod.method(
        "__cxxwrap_periodic",
        static_cast<std::tuple<bool,bool,bool> (*)(container_poly&)>(
            __cxxwrap_periodic
        )
    );
    mod.method(
        "__cxxwrap_point_inside",
        static_cast<bool (*)(container_poly&, double, double, double)>(
            __cxxwrap_point_inside
        )
    );
    mod.method(
        "__cxxwrap_draw_domain_gnuplot",
        static_cast<void (*)(container_poly&, FILE*)>(
            __cxxwrap_draw_domain_gnuplot_file
        )
    );
    mod.method(
        "__cxxwrap_draw_domain_gnuplot",
        static_cast<void (*)(container_poly&, const char*)>(
            __cxxwrap_draw_domain_gnuplot_filename
        )
    );
    mod.method(
        "__cxxwrap_draw_domain_pov",
        static_cast<void (*)(container_poly&, FILE*)>(
            __cxxwrap_draw_domain_pov_file
        )
    );
    mod.method(
        "__cxxwrap_draw_domain_pov",
        static_cast<void (*)(container_poly&, const char*)>(
            __cxxwrap_draw_domain_pov_filename
        )
    );
    mod.method(
        "__cxxwrap_total_particles",
        static_cast<int (*)(container_poly&)>(
            __cxxwrap_total_particles
        )
    );
    mod.method(
        "__cxxwrap_clear!",
        static_cast<void (*)(container_poly&)>(
            __cxxwrap_clear
        )
    );
    mod.method(
        "__cxxwrap_put!",
        static_cast<void (*)(container_poly&, int, double, double, double, double)>(
            [](container_poly& con, int i, double x, double y, double z, double r)
            {
                con.put(i, x, y, z, r);
            }
        )
    );
    mod.method(
        "__cxxwrap_put!",
        static_cast<void (*)(container_poly&, particle_order&, int, double, double, double, double)>(
            [](container_poly& con, particle_order& ord, int i, double x, double y, double z, double r)
            {
                con.put(ord, i, x, y, z, r);
            }
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container_poly&, FILE*)>(
            __cxxwrap_import_file
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container_poly&, particle_order&, FILE*)>(
            __cxxwrap_ordered_import_file
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container_poly&, const char*)>(
            __cxxwrap_import_filename
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container_poly&, particle_order&, const char*)>(
            __cxxwrap_ordered_import_filename
        )
    );
    mod.method(
        "__cxxwrap_compute_all_cells",
        static_cast<void (*)(container_poly&)>(
            __cxxwrap_compute_all_cells
        )
    );
    mod.method(
        "__cxxwrap_sum_cell_volumes",
        static_cast<void (*)(container_poly&)>(
            __cxxwrap_sum_cell_volumes
        )
    );
    mod.method(
        "__cxxwrap_draw_particles",
        static_cast<void (*)(container_poly&, FILE*)>(
            __cxxwrap_draw_particles_file
        )
    );
    mod.method(
        "__cxxwrap_draw_particles",
        static_cast<void (*)(container_poly&, const char*)>(
            __cxxwrap_draw_particles_filename
        )
    );
    mod.method(
        "__cxxwrap_draw_particles_pov",
        static_cast<void (*)(container_poly&, FILE*)>(
            __cxxwrap_draw_particles_pov_file
        )
    );
    mod.method(
        "__cxxwrap_draw_particles_pov",
        static_cast<void (*)(container_poly&, const char*)>(
            __cxxwrap_draw_particles_pov_filename
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_gnuplot",
        static_cast<void (*)(container_poly&, FILE*)>(
            __cxxwrap_draw_cells_file
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_gnuplot",
        static_cast<void (*)(container_poly&, const char*)>(
            __cxxwrap_draw_cells_filename
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_pov",
        static_cast<void (*)(container_poly&, FILE*)>(
            __cxxwrap_draw_cells_pov_file
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_pov",
        static_cast<void (*)(container_poly&, const char*)>(
            __cxxwrap_draw_cells_pov_filename
        )
    );
    mod.method(
        "__cxxwrap_find_voronoi_cell",
        static_cast<std::tuple<double, double, double, int, bool> (*)(container_poly&, double, double, double)>(
            __cxxwrap_find_voronoi_cell
        )
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container_poly&, c_loop_all&)>(
            __cxxwrap_compute_cell
        )
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container_poly&, c_loop_order&)>(
            __cxxwrap_compute_cell
        )
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container_poly&, c_loop_subset&)>(
            __cxxwrap_compute_cell
        )
    );
    mod.method(
        "__cxxwrap_point_inside_walls",
        static_cast<bool (*)(container_poly&, double, double, double)>(
            __cxxwrap_point_inside_walls
        )
    );
    mod.method(
        "__cxxwrap_add_wall!",
        static_cast<void (*)(container_poly&, wall_plane&)>(
            __cxxwrap_add_wall
        )
    );
    mod.method(
        "__cxxwrap_add_wall!",
        static_cast<void (*)(container_poly&, wall_sphere&)>(
            __cxxwrap_add_wall
        )
    );
    mod.method(
        "__cxxwrap_add_wall!",
        static_cast<void (*)(container_poly&, wall_cylinder&)>(
            __cxxwrap_add_wall
        )
    );
    mod.method(
        "__cxxwrap_add_wall!",
        static_cast<void (*)(container_poly&, wall_cone&)>(
            __cxxwrap_add_wall
        )
    );
    mod.method(
        "__cxxwrap_apply_walls!",
        static_cast<bool (*)(voronoicell_neighbor&, container_poly&, double, double, double)>(
            __cxxwrap_apply_walls
        )
    );

    // lambdas for walls
    mod.method(
        "__cxxwrap_point_inside",
        static_cast<bool (*)(wall_sphere&, double, double, double)>(__cxxwrap_point_inside)
    );

    mod.method(
        "__cxxwrap_point_inside",
        static_cast<bool (*)(wall_cylinder&, double, double, double)>(__cxxwrap_point_inside)
    );

    mod.method(
        "__cxxwrap_point_inside",
        static_cast<bool (*)(wall_cone&, double, double, double)>(__cxxwrap_point_inside)
    );

    mod.method(
        "__cxxwrap_point_inside",
        static_cast<bool (*)(wall_plane&, double, double, double)>(__cxxwrap_point_inside)
    );

    // Class Container Periodic

    mod.method(
        "__cxxwrap_bounds",
        static_cast<std::tuple<double,double,double,double,double,double> (*)(container_periodic&)>(
            __cxxwrap_bounds_triclinic
        )
    );
    mod.method(
        "__cxxwrap_clear!",
        static_cast<void (*)(container_periodic&)>(
            __cxxwrap_clear
        )
    );
    mod.method(
        "__cxxwrap_put!",
        static_cast<void (*)(container_periodic&, int, double, double, double)>(
            [](container_periodic& con, int i, double x, double y, double z)
            {
                con.put(i, x, y, z);
            }
        )
    );
    mod.method(
        "__cxxwrap_put!",
        static_cast<void (*)(container_periodic&, particle_order&, int, double, double, double)>(
            [](container_periodic& con, particle_order& ord, int i, double x, double y, double z)
            {
                con.put(ord, i, x, y, z);
            }
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container_periodic&, FILE*)>(
            __cxxwrap_import_file
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container_periodic&, particle_order&, FILE*)>(
            __cxxwrap_ordered_import_file
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container_periodic&, const char*)>(
            __cxxwrap_import_filename
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container_periodic&, particle_order&, const char*)>(
            __cxxwrap_ordered_import_filename
        )
    );
    mod.method(
        "__cxxwrap_compute_all_cells",
        static_cast<void (*)(container_periodic&)>(
            __cxxwrap_compute_all_cells
        )
    );
    mod.method(
        "__cxxwrap_sum_cell_volumes",
        static_cast<void (*)(container_periodic&)>(
            __cxxwrap_sum_cell_volumes
        )
    );
    mod.method(
        "__cxxwrap_draw_particles",
        static_cast<void (*)(container_periodic&, FILE*)>(
            __cxxwrap_draw_particles_file
        )
    );
    mod.method(
        "__cxxwrap_draw_particles",
        static_cast<void (*)(container_periodic&, const char*)>(
            __cxxwrap_draw_particles_filename
        )
    );
    mod.method(
        "__cxxwrap_draw_particles_pov",
        static_cast<void (*)(container_periodic&, FILE*)>(
            __cxxwrap_draw_particles_pov_file
        )
    );
    mod.method(
        "__cxxwrap_draw_particles_pov",
        static_cast<void (*)(container_periodic&, const char*)>(
            __cxxwrap_draw_particles_pov_filename
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_gnuplot",
        static_cast<void (*)(container_periodic&, FILE*)>(
            __cxxwrap_draw_cells_file
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_gnuplot",
        static_cast<void (*)(container_periodic&, const char*)>(
            __cxxwrap_draw_cells_filename
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_pov",
        static_cast<void (*)(container_periodic&, FILE*)>(
            __cxxwrap_draw_cells_pov_file
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_pov",
        static_cast<void (*)(container_periodic&, const char*)>(
            __cxxwrap_draw_cells_pov_filename
        )
    );
    mod.method(
        "__cxxwrap_find_voronoi_cell",
        static_cast<std::tuple<double, double, double, int, bool> (*)(container_periodic&, double, double, double)>(
            __cxxwrap_find_voronoi_cell
        )
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container_periodic&, c_loop_all_periodic&)>(
            __cxxwrap_compute_cell
        )
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container_periodic&, c_loop_order_periodic&)>(
            __cxxwrap_compute_cell
        )
    );

    // Class Container Periodic Poly

    mod.method(
        "__cxxwrap_bounds",
        static_cast<std::tuple<double,double,double,double,double,double> (*)(container_periodic_poly&)>(
            __cxxwrap_bounds_triclinic
        )
    );
    mod.method(
        "__cxxwrap_clear!",
        static_cast<void (*)(container_periodic_poly&)>(
            __cxxwrap_clear
        )
    );
    mod.method(
        "__cxxwrap_put!",
        static_cast<void (*)(container_periodic_poly&, int, double, double, double, double)>(
            [](container_periodic_poly& con, int i, double x, double y, double z, double r)
            {
                con.put(i, x, y, z, r);
            }
        )
    );
    mod.method(
        "__cxxwrap_put!",
        static_cast<void (*)(container_periodic_poly&, particle_order&, int, double, double, double, double)>(
            [](container_periodic_poly& con, particle_order& ord, int i, double x, double y, double z, double r)
            {
                con.put(ord, i, x, y, z, r);
            }
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container_periodic_poly&, FILE*)>(
            __cxxwrap_import_file
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container_periodic_poly&, particle_order&, FILE*)>(
            __cxxwrap_ordered_import_file
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container_periodic_poly&, const char*)>(
            __cxxwrap_import_filename
        )
    );
    mod.method(
        "__cxxwrap_import!",
        static_cast<void (*)(container_periodic_poly&, particle_order&, const char*)>(
            __cxxwrap_ordered_import_filename
        )
    );
    mod.method(
        "__cxxwrap_compute_all_cells",
        static_cast<void (*)(container_periodic_poly&)>(
            __cxxwrap_compute_all_cells
        )
    );
    mod.method(
        "__cxxwrap_sum_cell_volumes",
        static_cast<void (*)(container_periodic_poly&)>(
            __cxxwrap_sum_cell_volumes
        )
    );
    mod.method(
        "__cxxwrap_draw_particles",
        static_cast<void (*)(container_periodic_poly&, FILE*)>(
            __cxxwrap_draw_particles_file
        )
    );
    mod.method(
        "__cxxwrap_draw_particles",
        static_cast<void (*)(container_periodic_poly&, const char*)>(
            __cxxwrap_draw_particles_filename
        )
    );
    mod.method(
        "__cxxwrap_draw_particles_pov",
        static_cast<void (*)(container_periodic_poly&, FILE*)>(
            __cxxwrap_draw_particles_pov_file
        )
    );
    mod.method(
        "__cxxwrap_draw_particles_pov",
        static_cast<void (*)(container_periodic_poly&, const char*)>(
            __cxxwrap_draw_particles_pov_filename
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_gnuplot",
        static_cast<void (*)(container_periodic_poly&, FILE*)>(
            __cxxwrap_draw_cells_file
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_gnuplot",
        static_cast<void (*)(container_periodic_poly&, const char*)>(
            __cxxwrap_draw_cells_filename
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_pov",
        static_cast<void (*)(container_periodic_poly&, FILE*)>(
            __cxxwrap_draw_cells_pov_file
        )
    );
    mod.method(
        "__cxxwrap_draw_cells_pov",
        static_cast<void (*)(container_periodic_poly&, const char*)>(
            __cxxwrap_draw_cells_pov_filename
        )
    );
    mod.method(
        "__cxxwrap_find_voronoi_cell",
        static_cast<std::tuple<double, double, double, int, bool> (*)(container_periodic_poly&, double, double, double)>(
            __cxxwrap_find_voronoi_cell
        )
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container_periodic_poly&, c_loop_all_periodic&)>(
            __cxxwrap_compute_cell
        )
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container_periodic_poly&, c_loop_order_periodic&)>(
            __cxxwrap_compute_cell
        )
    );

    mod.method(
        "__cxxwrap_compute_ghost_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container&, double, double, double)>(
            __cxxwrap_compute_ghost_cell
        )
    );
    mod.method(
        "__cxxwrap_compute_ghost_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container_poly&, double, double, double, double)>(
            __cxxwrap_compute_ghost_cell_poly
        )
    );
    mod.method(
        "__cxxwrap_compute_ghost_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container_periodic&, double, double, double)>(
            __cxxwrap_compute_ghost_cell
        )
    );
    mod.method(
        "__cxxwrap_compute_ghost_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container_periodic_poly&, double, double, double, double)>(
            __cxxwrap_compute_ghost_cell_poly
        )
    );

}
