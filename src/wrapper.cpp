#include "jlcxx/jlcxx.hpp"
#include "jlcxx/tuple.hpp"
#include "jlcxx/stl.hpp"
#include "voro++.hh"

/* # General notes
   Exported names starting with __cxxwrap_ mean that the return type or signature 
   is to be altered in the exported Julia code.
   ## Examples
   
   * The Julia function will return the modified cell:

   void compute_cell(voronoicell& vc, container& con, c_loop_all& cl) {
        con.compute_cell(vc, cl);
   }

   * bool is wrapped as CxxBool, needs conversion to Bool to use in logical expressions

   bool is_periodic_x(container& con) {
        return con.xperiodic;
   }
*/

/* Get and Set methods for accessing public members from container class container.hh */
// int get_particle_id(voro::container& con, int i, int j) {
//     return con.id[i][j];
// }

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    using namespace voro;

    // Class Container
    auto class_container = mod.add_type<container>(
        "RawContainer", jlcxx::julia_type("AbstractRawContainer", "VoroPlusPlus")
    )
        .constructor<double, double, double, double, double, double, int, int, int, bool, bool, bool, int>()
    ;

    // Class Container Poly
    auto class_container_poly = mod.add_type<container_poly>(
        "RawContainerPoly", jlcxx::julia_type("AbstractRawContainer", "VoroPlusPlus")
    )
        .constructor<double, double, double, double, double, double, int, int, int, bool, bool, bool, int>()
    ;

    // Class Containter Periodic Poly
    auto class_container_periodic_poly = mod.add_type<container_periodic_poly>("RawContainerTriclinic")
        .constructor<double, double, double, double, double, double, int, int, int, int>()
    ;

    // Class VoronoiCell
    auto class_voronoicell_neighbor = mod.add_type<voronoicell_neighbor>(
        "VoronoiCell", jlcxx::julia_type("AbstractVoronoiCell", "VoroPlusPlus")
    )
        .constructor<>()
        .constructor<double>()
        .constructor<container&>()
        .constructor<container_poly&>()
    ;

    // Class Particle Order
    auto class_particle_order = mod.add_type<particle_order>(
        "InsertionOrder",  jlcxx::julia_type("ContainerIterationOrder", "VoroPlusPlus")
    )
        .constructor<>()
    ;

    // Class Container Iterator
    auto class_c_loop_all = mod.add_type<c_loop_all>("RawContainerIterator")
        .constructor<container&>()
        .constructor<container_poly&>()
    ;

    // Class Container Iterator Subset
    auto class_c_loop_subset = mod.add_type<c_loop_subset>("ContainerSubsetIterator")
        .constructor<container&>()
        .constructor<container_poly&>()
    ;

    // Class Container Iterator Order
    auto class_c_loop_order = mod.add_type<c_loop_order>("InsertionOrderIterator")
        .constructor<container&, particle_order&>()
        .constructor<container_poly&, particle_order&>()
    ;

    // Class Wall   
    auto class_wall_sphere = mod.add_type<wall_sphere>(
        "WallSphere", jlcxx::julia_type("AbstractWall", "VoroPlusPlus")
    )
        .constructor<double, double, double, double>()
        .constructor<double, double, double, double, int>()
    ;

    auto class_wall_cylinder = mod.add_type<wall_cylinder>(
        "WallCylinder", jlcxx::julia_type("AbstractWall", "VoroPlusPlus")
    )
        .constructor<double, double, double, double, double, double, double, int>()
        .constructor<double, double, double, double, double, double, double>()
    ;

    auto class_wall_cone = mod.add_type<wall_cone>(
        "WallCone", jlcxx::julia_type("AbstractWall", "VoroPlusPlus")
    )
        .constructor<double, double, double, double, double, double, double, int>()
        .constructor<double, double, double, double, double, double, double>()
    ;

    auto class_wall_plane = mod.add_type<wall_plane>(
        "WallPlane", jlcxx::julia_type("AbstractWall", "VoroPlusPlus")
    )
        .constructor<double, double, double, double, int>()
        .constructor<double, double, double, double>()
    ;

    // container_base methods
    
    auto __cxxwrap_bounds = [] (auto &con)
    {
        return std::make_tuple(con.ax, con.ay, con.az, con.bx, con.by, con.bz);
    };

    auto __cxxwrap_periodic = [] (auto &con)
    {
        return std::make_tuple(con.xperiodic, con.yperiodic, con.zperiodic);
    };

    auto __cxxwrap_point_inside = [] (auto &w, double x, double y, double z) 
    {
        return w.point_inside(x, y, z);
    };

    auto __cxxwrap_point_inside_walls = [] (auto &con, double x, double y, double z) 
    {
        return con.point_inside_walls(x, y, z);
    };

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

    auto __cxxwrap_ordered_import_filename = [] (auto& con, particle_order& ord, const char* fp)
    {
        con.import(ord, filename);
    };
    
    auto __cxwrap_compute_all_cells = [] (auto& con) {con.compute_all_cells()};

    auto __cxwrap_sum_cell_volumes = [] (auto& con) {con.sum_cell_volumes()};
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
        "__cxxwrap_point_inside_walls",
        static_cast<bool (*)(container&, double, double, double)>(
            __cxxwrap_point_inside_walls
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
        static_cast<void (*)(container&, int, double, double, double)>([](container& con, int i, double x, double y, double z)
        {
            con.put(i, x, y, z);
        })
    );
    mod.method(
        "__cxxwrap_put!",
        static_cast<void (*)(container&, particle_order&, int, double, double, double)>([](container& con, particle_order& ord, int i, double x, double y, double z)
        {
            con.put(ord, i, x, y, z);
        })
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
        static_cast<double (*)(container&)>(
            __cxxwrap_sum_cell_volumes
        )
    );
    mod.method(
        "__cxxwrap_draw_particles",
        static_cast<double (*)(container&, FILE*)>(
            __cxxwrap_draw_particles_file
        )
    );
    mod.method(
        "__cxxwrap_draw_particles",
        static_cast<double (*)(container&, const char*)>(
            __cxxwrap_draw_particles_filename
        )
    );
    mod.method(
        "__cxxwrap_draw_particles_pov",
        static_cast<double (*)(container&, FILE*)>(
            __cxxwrap_draw_particles_pov_file
        )
    );
    mod.method(
        "__cxxwrap_draw_particles_pov",
        static_cast<double (*)(container&, const char*)>(
            __cxxwrap_draw_particles_pov_filename
        )
    );

    
    class_container
        .method(
            "draw_particles_pov",
            static_cast<void (container::*)(const char*)>(&container::draw_particles_pov)
        )
        .method(
            "draw_particles_pov",
            static_cast<void (container::*)(FILE*)>(&container::draw_particles_pov)
        )
        .method(
            "draw_cells_gnuplot",
            static_cast<void (container::*)(const char*)>(&container::draw_cells_gnuplot)
        )
        .method(
            "draw_cells_gnuplot",
            static_cast<void (container::*)(FILE*)>(&container::draw_cells_gnuplot)
        )
        .method(
            "draw_cells_pov",
            static_cast<void (container::*)(const char*)>(&container::draw_cells_pov)
        )
        .method(
            "draw_cells_pov",
            static_cast<void (container::*)(FILE*)>(&container::draw_cells_pov)
        )
        .method("draw_domain_gnuplot",
            static_cast<void (container::*)(const char*)>(&container::draw_domain_gnuplot)
        )
        .method("draw_domain_pov",
            static_cast<void (container::*)(const char*)>(&container::draw_domain_pov)
        )
        ;


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

    class_container_poly
        .method(
            "__cxxwrap_put!",
            static_cast<void (container_poly::*)(int, double, double, double, double)>(&container_poly::put)
        )
        .method(
            "__cxxwrap_put!",
            static_cast<void (container_poly::*)(particle_order&, int, double, double, double, double)>(&container_poly::put)
        )
        .method(
            "__cxxwrap_clear!",
            &container_poly::clear
        )
        .method(
            "__cxxwrap_isinside", 
            static_cast<bool (container_poly::*)(double, double, double)>(&container_poly::point_inside)
        )
        .method(
            "compute_all_cells",
            &container_poly::compute_all_cells
        )
        .method(
            "sum_cell_volumes", 
            &container_poly::sum_cell_volumes
        )
        .method("import!", static_cast<void (container_poly::*)(const char*)>(&container_poly::import))
        .method("draw_particles", static_cast<void (container_poly::*)(const char*)>(&container_poly::draw_particles))
        .method("draw_particles_pov", static_cast<void (container_poly::*)(const char*)>(&container_poly::draw_particles_pov))
        .method("draw_cells_gnuplot", static_cast<void (container_poly::*)(const char*)>(&container_poly::draw_cells_gnuplot))
        .method("draw_cells_pov", static_cast<void (container_poly::*)(const char*)>(&container_poly::draw_cells_pov))
        .method("import!", static_cast<void (container_poly::*)(FILE*)>(&container_poly::import))
        .method("draw_particles", static_cast<void (container_poly::*)(FILE*)>(&container_poly::draw_particles))
        .method("draw_particles_pov", static_cast<void (container_poly::*)(FILE*)>(&container_poly::draw_particles_pov))
        .method("draw_cells_gnuplot", static_cast<void (container_poly::*)(FILE*)>(&container_poly::draw_cells_gnuplot))
        .method("draw_cells_pov", static_cast<void (container_poly::*)(FILE*)>(&container_poly::draw_cells_pov))
        ;

    class_c_loop_all
        .method(
            "__cxxwrap_start!",
            &c_loop_all::start
        )
        .method(
            "__cxxwrap_inc!",
            &c_loop_all::inc
        )
        ;
    
    class_c_loop_subset
        .method(
            "__cxxwrap_start!",
            &c_loop_subset::start
        )
        .method(
            "__cxxwrap_inc!",
            &c_loop_subset::inc
        )
        .method(
            "__cxxwrap_set_bounds_to_sphere!",
            static_cast<void (c_loop_subset::*)(double, double, double, double, bool)>(&c_loop_subset::setup_sphere)
        )
        .method(
            "__cxxwrap_set_bounds_to_box!",
            static_cast<void (c_loop_subset::*)(double, double, double, double, double, double, bool)>(&c_loop_subset::setup_box)
        )
    ;


    class_c_loop_order
        .method(
            "__cxxwrap_start!",
            &c_loop_order::start
        )
        .method(
            "__cxxwrap_inc!",
            &c_loop_order::inc
        )
        ;


    // Class VoronoiCell
    class_voronoicell_neighbor
        .method(
            "__cxxwrap_init!",
            &voronoicell_neighbor::init
        )
        .method(
            "__cxxwrap_init_octahedron!",
            &voronoicell_neighbor::init_octahedron
        )
        .method(
            "__cxxwrap_init_tetrahedron!",
            &voronoicell_neighbor::init_tetrahedron
        )
        .method(
            "__cxxwrap_nplane!",
            static_cast<bool (voronoicell_neighbor::*)(double, double, double, double, int)>(&voronoicell_neighbor::nplane)
        )
        .method(
            "__cxxwrap_nplane!",
            static_cast<bool (voronoicell_neighbor::*)(double, double, double, int)>(&voronoicell_neighbor::nplane)
        )
        .method("volume", &voronoicell_neighbor::volume)
        //.method("check_relations", &voronoicell_neighbor::check_relations)
        //.method("check_duplicates", &voronoicell_neighbor::check_duplicates)
        .method("max_radius_squared", &voronoicell_neighbor::max_radius_squared)
        .method("number_of_edges", &voronoicell_neighbor::number_of_edges)
        .method("total_edge_distance", &voronoicell_neighbor::total_edge_distance)
        .method("number_of_faces", &voronoicell_neighbor::number_of_faces)
        .method("surface_area", &voronoicell_neighbor::surface_area)
        .method(
            "__cxxwrap_draw_gnuplot",
            static_cast<void (voronoicell_neighbor::*)(double, double, double, const char*)>(
                &voronoicell_neighbor::draw_gnuplot
            )
        )
        .method(
            "__cxxwrap_draw_gnuplot",
            static_cast<void (voronoicell_neighbor::*)(double, double, double, FILE*)>(
                &voronoicell_neighbor::draw_gnuplot
            )
        )
        .method("__cxxwrap_draw_pov", static_cast<void (voronoicell_neighbor::*)(double, double, double, const char*)>(&voronoicell_neighbor::draw_pov))
        .method("__cxxwrap_draw_pov", static_cast<void (voronoicell_neighbor::*)(double, double, double, FILE*)>(&voronoicell_neighbor::draw_pov))
        .method("__cxxwrap_draw_pov_mesh", static_cast<void (voronoicell_neighbor::*)(double, double, double, const char*)>(&voronoicell_neighbor::draw_pov_mesh))
        .method("__cxxwrap_draw_pov_mesh", static_cast<void (voronoicell_neighbor::*)(double, double, double, FILE*)>(&voronoicell_neighbor::draw_pov_mesh))
        .method("__cxxwrap_plane_intersects", static_cast<bool (voronoicell_neighbor::*)(double, double, double, double)>(&voronoicell_neighbor::plane_intersects))
        // Inherited from voronoicell_base
        .method("__cxxwrap_translate!", &voronoicell_neighbor::translate)
        .method("print_edges", &voronoicell_neighbor::print_edges)
        .method("__cycle_up", static_cast<int (voronoicell_neighbor::*)(int, int)>(&voronoicell_neighbor::cycle_up))
        .method("__cycle_down", static_cast<int (voronoicell_neighbor::*)(int, int)>(&voronoicell_neighbor::cycle_down))
        ;

    // Class Containter Periodic Poly (conprdply)
    class_container_periodic_poly
        .method("__cxxwrap_put!", static_cast<void (container_periodic_poly::*)(int, double, double, double, double)>(&container_periodic_poly::put))
        // Type mismatch from Voro++
        //.method("conprdply_compute_ghost_cell", static_cast<bool (container_periodic_poly::*)(voronoicell&, double, double, double, double)>(&container_periodic_poly::compute_ghost_cell))
        .method("__cxxwrap_draw_particles", static_cast<void (container_periodic_poly::*)(const char*)>(&container_periodic_poly::draw_particles))
        .method("__cxxwrap_draw_cells_gnuplot", static_cast<void (container_periodic_poly::*)(const char*)>(&container_periodic_poly::draw_cells_gnuplot))
        .method("__cxxwrap_draw_domain_gnuplot", static_cast<void (container_periodic_poly::*)(const char*)>(&container_periodic_poly::draw_domain_gnuplot))
        ;

    // lambdas for walls

    auto __cxxwrap_add_wall = [] (auto &con, auto &wall) { con.add_wall(wall); };
            // Type mismatch from Voro++
        //.method("contains_neighbor", static_cast<void (container::*)(const char*)>(&container::contains_neighbor))
    mod.method("__cxxwrap_add_wall!",
        static_cast<void (*)(container&, wall_plane&)>(__cxxwrap_add_wall)
    );
    mod.method("__cxxwrap_add_wall!",
        static_cast<void (*)(container&, wall_sphere&)>(__cxxwrap_add_wall)
    );
    mod.method("__cxxwrap_add_wall!",
        static_cast<void (*)(container&, wall_cylinder&)>(__cxxwrap_add_wall)
    );
    mod.method("__cxxwrap_add_wall!",
        static_cast<void (*)(container&, wall_cone&)>(__cxxwrap_add_wall)
    );
    mod.method("__cxxwrap_add_wall!",
        static_cast<void (*)(container_poly&, wall_plane&)>(__cxxwrap_add_wall)
    );
    mod.method("__cxxwrap_add_wall!",
        static_cast<void (*)(container_poly&, wall_sphere&)>(__cxxwrap_add_wall)
    );
    mod.method("__cxxwrap_add_wall!",
        static_cast<void (*)(container_poly&, wall_cylinder&)>(__cxxwrap_add_wall)
    );
    mod.method("__cxxwrap_add_wall!",
        static_cast<void (*)(container_poly&, wall_cone&)>(__cxxwrap_add_wall)
    );

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

    mod.method(
        "__cxxwrap_point_inside_walls",
        static_cast<bool (*)(container&, double, double, double)>(__cxxwrap_point_inside_walls)
    );
    mod.method(
        "__cxxwrap_point_inside_walls",
        static_cast<bool (*)(container_poly&, double, double, double)>(__cxxwrap_point_inside_walls)
    );

    // mod.method(
    //     "__cxxwrap_point_inside",
    //     static_cast<bool (*)(double, double, double, wall&)>(
    //         [] (double x, double y, double z, wall& w) {
    //             return w.point_inside(x, y, z);
    //         }
    //     )
    // );

    // mod.method(
    //     "__cxxwrap_cut_cell!",
    //     static_cast<bool (*)(voronoicell_neighbor&, double, double, double, wall&)>(
    //         [] (voronoicell_neighbor& vc, double x, double y, double z, wall& w) {
    //             return w.cut_cell(vc, x, y, z);
    //         }
    //     )
    // );

    // lambdas for container and container_poly

    auto __cxxwrap_find_cell = [] (auto &con, double x, double y, double z)
    {
        int pid;
        double rx, ry, rz;
        bool found = con.find_voronoi_cell(x, y, z, rx, ry, rz, pid);
        return std::make_tuple(found, pid, rx, ry, rz);
    };
    
    mod.method(
        "__cxxwrap_find_cell",
        static_cast<std::tuple<bool,int,double,double,double> (*)(container&, double, double, double)>(
            __cxxwrap_find_cell
        )
    );
    mod.method(
        "__cxxwrap_find_cell",
        static_cast<std::tuple<bool,int,double,double,double> (*)(container_poly&, double, double, double)>(
            __cxxwrap_find_cell
        )
    );

    // lambdas for voronoicell
    mod.method(
        "__cxxwrap_centroid",
        static_cast<std::tuple<double,double,double> (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) {
                double x, y, z;
                v.centroid(x, y, z);
                return std::make_tuple(x, y, z);
            }
        )
    );

    // public attribute accessors
    mod.method(
        "__get_tol",
        static_cast<double (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.tol; }
        )
    );
    mod.method(
        "__get_current_vertices",
        static_cast<int (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.current_vertices; }
        )
    );
    mod.method(
        "__get_p",
        static_cast<int (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.p; }
        )
    );
    mod.method(
        "__get_up",
        static_cast<int (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.up; }
        )
    );
    mod.method(
        "__cxxwrap_get_ed",
        static_cast<int** (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.ed; }
        )
    );
    mod.method(
        "__cxxwrap_get_nu",
        static_cast<int* (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.nu; }
        )
    );
    mod.method(
        "__cxxwrap_get_pts",
        static_cast<double* (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.pts; }
        )
    );
    mod.method(
        "__cxxwrap_get_mne",
        static_cast<int** (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.mne; }
        )
    );
    mod.method(
        "__cxxwrap_get_ne",
        static_cast<int** (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.ne; }
        )
    );
    mod.method(
        "__get_ed_ij",
        static_cast<int (*)(voronoicell_neighbor&, size_t, size_t)>(
            [] (voronoicell_neighbor& v, size_t i, size_t j) { return v.ed[i][j]; }
        )
    );
    mod.method(
        "__set_ed_ij!",
        static_cast<int (*)(voronoicell_neighbor&, size_t, size_t, int)>(
            [] (voronoicell_neighbor& v, size_t i, size_t j, int k) { v.ed[i][j] = k; return k; }
        )
    );
    mod.method("__cxxwrap_copy!",
        static_cast<void (*)(voronoicell_neighbor&, voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& dest, voronoicell_neighbor& src) {
                dest = src;
            }
        )
    );
    mod.method(
        "__cxxwrap_get_vertex_orders!",
        static_cast<void (*)(std::vector<int>&, voronoicell_neighbor&)>(
            [] (std::vector<int>& v, voronoicell_neighbor& vc) { vc.vertex_orders(v); }   
        )
    );
    mod.method(
        "__cxxwrap_get_neighbors!",
        static_cast<void (*)(std::vector<int>&, voronoicell_neighbor&)>(
            [] (std::vector<int>& v, voronoicell_neighbor& vc) { vc.neighbors(v); }   
        )
    );
    mod.method(
        "__cxxwrap_vertices!",
        static_cast<void (*)(std::vector<double>&, voronoicell_neighbor&)>(
            [] (std::vector<double>& v, voronoicell_neighbor& vc) { vc.vertices(v); }   
        )
    );
    mod.method(
        "__cxxwrap_vertices!",
        static_cast<void (*)(std::vector<double>&, voronoicell_neighbor&, double, double, double)>(
            [] (std::vector<double>& v, voronoicell_neighbor& vc, double x, double y, double z)
            {
                vc.vertices(x, y, z, v);
            }   
        )
    );
    mod.method(
        "__cxxwrap_face_perimeters!",
        static_cast<void (*)(std::vector<double>&, voronoicell_neighbor&)>(
            [](std::vector<double>& v, voronoicell_neighbor& vc)
            {
                vc.face_perimeters(v);
            }
        )
    );
    mod.method(
        "__cxxwrap_face_areas!",
        static_cast<void (*)(std::vector<double>&, voronoicell_neighbor&)>(
            [](std::vector<double>& v, voronoicell_neighbor& vc)
            {
                vc.face_areas(v);
            }
        )
    );
    mod.method(
        "__cxxwrap_normals!",
        static_cast<void (*)(std::vector<double>&, voronoicell_neighbor&)>(
            [](std::vector<double>& v, voronoicell_neighbor& vc)
            {
                vc.normals(v);
            }
        )
    );
    mod.method(
        "__cxxwrap_face_vertices!",
        static_cast<void (*)(std::vector<int>&, voronoicell_neighbor&)>(
            [](std::vector<int>& v, voronoicell_neighbor& vc)
            {
                vc.face_vertices(v);
            }
        )
    );
    mod.method(
        "__cxxwrap_face_orders!",
        static_cast<void (*)(std::vector<int>&, voronoicell_neighbor&)>(
            [](std::vector<int>& v, voronoicell_neighbor& vc)
            {
                vc.face_orders(v);
            }
        )
    );
    mod.method(
        "__cxxwrap_face_freq_table!",
        static_cast<void (*)(std::vector<int>&, voronoicell_neighbor&)>(
            [](std::vector<int>& v, voronoicell_neighbor& vc)
            {
                vc.face_freq_table(v);
            }
        )
    );

    // lambdas for loops

    auto pos = [] (auto &cl) 
    {
        double x, y, z;
        cl.pos(x, y, z);
        return std::make_tuple(x, y, z);
    };

    auto particle_info = [] (auto &cl)
    {
        int pid;
        double x, y, z, r;
        cl.pos(pid, x, y, z, r);
        return std::make_tuple(pid, x, y, z, r);
    };

    mod.method(
        "pos",
        static_cast<std::tuple<double,double,double> (*)(c_loop_all&)>(pos)
    );
    mod.method(
        "pos",
        static_cast<std::tuple<double,double,double> (*)(c_loop_order&)>(pos)
    );
    mod.method(
        "pos",
        static_cast<std::tuple<double,double,double> (*)(c_loop_subset&)>(pos)
    );
    mod.method(
        "particle_info",
        static_cast<std::tuple<int,double,double,double,double> (*)(c_loop_all&)>(particle_info)
    );
    mod.method(
        "particle_info",
        static_cast<std::tuple<int,double,double,double,double> (*)(c_loop_order&)>(particle_info)
    );
    mod.method(
        "particle_info",
        static_cast<std::tuple<int,double,double,double,double> (*)(c_loop_subset&)>(particle_info)
    );

    // Anonymus functions for special cases when two types are needed
    auto compute_cell = [] (voronoicell_neighbor& vc, auto& con, auto& itr)
    {
        return con.compute_cell(vc, itr);
    };
    
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container&, c_loop_all&)>(compute_cell)
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container&, c_loop_order&)>(compute_cell)
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container_poly&, c_loop_all&)>(compute_cell)
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container_poly&, c_loop_order&)>(compute_cell)
    );

    mod.method("apply_walls!", [] (voronoicell_neighbor& vc, container& con, double x, double y, double z){ return con.apply_walls(vc, x, y, z);});

    // Public Menbers from Container Class
    // mod.method("get_particle_id", &get_particle_id);
}