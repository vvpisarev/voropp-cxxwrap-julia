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
int get_particle_id(voro::container& con, int i, int j) {

    return con.id[i][j];
}



/* Extern Method for special cases Voronoicell */
extern "C" {

    void draw_gnuplot_voronoicell(voro::voronoicell_neighbor* vc, double x, double y, double z, FILE *fp) { vc->draw_gnuplot(x, y, z, fp); }

    // Version implemented from Voro++  >>>  void output_vertices(FILE *fp=stdout), ignoring stdout by default
    // without position parameters
    void output_vertices_nopos(voro::voronoicell_neighbor* vc, FILE *fp) { vc->output_vertices(fp); }

    // Version implemented from Voro++  >>>  void output_vertices(double x,double y,double z,FILE *fp=stdout), ignoring stdout by default
    // with position parameters
    void output_vertices_positions(voro::voronoicell_neighbor* vc, double x, double y, double z, FILE *fp) { vc->output_vertices(x, y, z, fp); }

    // Version implemented from Voro++  >>>  void output_vertex_orders(FILE *fp=stdout), ignoring stdout by default
    void output_vertex_orders_vorocell(voro::voronoicell_neighbor* vc, FILE *fp) { vc->output_vertex_orders(fp); }

    // Version implemented from Voro++  >>>  inline void output_face_perimeters(FILE *fp=stdout), ignoring stdout by default
    void output_face_perimeters_vorocell(voro::voronoicell_neighbor* vc, FILE *fp) { vc->output_face_perimeters(fp); }

    // Version implemented from Voro++  >>>  inline void output_face_freq_table(FILE *fp=stdout), ignoring stdout by default
    void output_face_freq_table_vorocell(voro::voronoicell_neighbor* vc, FILE *fp) { vc->output_face_freq_table(fp); }

    // Version implemented from Voro++  >>>  inline void output_face_orders(FILE *fp=stdout), ignoring stdout by default
    void output_face_orders_vorocell(voro::voronoicell_neighbor* vc, FILE *fp) { vc->output_face_orders(fp); }
    
    // Version implemented from Voro++  >>>  inline void output_face_areas(FILE *fp=stdout), ignoring stdout by default
    void output_face_areas_vorocell(voro::voronoicell_neighbor* vc, FILE *fp) { vc->output_face_areas(fp); }
    
    // Version implemented from Voro++  >>>  inline void output_normals(FILE *fp=stdout), ignoring stdout by default
    void output_normals_vorocell(voro::voronoicell_neighbor* vc, FILE *fp) { vc->output_normals(fp); }
    
    // Version implemented from Voro++  >>>  inline void output_face_vertices(FILE *fp=stdout), ignoring stdout by default
    void output_face_vertices_vorocell(voro::voronoicell_neighbor* vc, FILE *fp) { vc->output_face_vertices(fp); }

    bool compute_ghost_cell_conprdply(voro::container_periodic_poly* con, voro::voronoicell_neighbor* c, double x, double y, double z, double r) {

        return con->compute_ghost_cell(c, x, y, z, r);
    }

}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    using namespace voro;

    // Class Particle Order
    mod.add_type<particle_order>(
        "InsertionOrder",  jlcxx::julia_type("ContainerIterationOrder", "VoroPlusPlus")
    )
        .constructor<>()
    ;

    // Class Wall   
    mod.add_type<wall>("Wall")
        .method("point_inside", static_cast<bool (wall::*)(double,double,double)>(&wall::point_inside))
        //.method("cut_cell", &wall::cut_cell)
        ;
    
    
    // Class Wall_List
    mod.add_type<wall_list>("Wall_List")
        .constructor<>()
        .method("add_wall", static_cast<void (wall_list::*)(wall&)>(&wall_list::add_wall))
        ;


    // Class Container
    mod.add_type<container>(
        "RawContainer", jlcxx::julia_type("AbstractRawContainer", "VoroPlusPlus")
    )
        .constructor<double, double, double, double, double, double, int, int, int, bool, bool, bool, int>()
        .method(
            "__cxxwrap_put!",
            static_cast<void (container::*)(int, double, double, double)>(&container::put)
        )
        .method(
            "__cxxwrap_put!",
            static_cast<void (container::*)(particle_order&, int, double, double, double)>(&container::put)
        )
        .method(
            "import!",
            static_cast<void (container::*)(const char*)>(&container::import)
        )
        // When exporting to Julia, it is no possible to derive automatically type FILE
        //.method("draw_cells_gnuplot", static_cast<void (container::*)(FILE*)>(&container::draw_cells_gnuplot))
        .method(
            "draw_particles",
            static_cast<void (container::*)(const char*)>(&container::draw_particles)
        )
        .method(
            "draw_particles_pov",
            static_cast<void (container::*)(const char*)>(&container::draw_particles_pov)
        )
        .method(
            "draw_cells_gnuplot",
            static_cast<void (container::*)(const char*)>(&container::draw_cells_gnuplot)
        )
        .method(
            "print_custom",
            static_cast<void (container::*)(const char*, const char*)>(&container::print_custom)
        )
        .method(
            "draw_cells_pov",
            static_cast<void (container::*)(const char*)>(&container::draw_cells_pov)
        )
        .method(
            "__cxxwrap_clear!",
            &container::clear
        )
        .method(
            "compute_all_cells",
            &container::compute_all_cells
        )
        .method(
            "sum_cell_volumes", 
            &container::sum_cell_volumes
        )
        .method(
            "__cxxwrap_isinside", 
            static_cast<bool (container::*)(double, double, double)>(&container::point_inside)
        )
        .method(
            "region_count",
            &container::region_count
        )
        // When exporting to Julia, it is needed to defien first voronoicell type
        //.method("initialize_voronoicell", static_cast<bool (container::*)(voro::voronoicell&, int, int, int, int, int, int&, int&, int&, double&, double&, double&, int&)>(&container::initialize_voronoicell))
        .method(
            "initialize_search",
            static_cast<void (container::*)(int, int, int, int, int&, int&, int&, int&)>(&container::initialize_search)
        )
        .method(
            "frac_pos",
            static_cast<void (container::*)(double, double, double, double, double, double, double&, double&, double&)>(&container::frac_pos)
        )
        .method(
            "region_index",
            static_cast<int (container::*)(int, int, int, int, int, int, double&, double&, double&, int&)>(&container::region_index)
        )
        .method("draw_domain_gnuplot",
            static_cast<void (container::*)(const char*)>(&container::draw_domain_gnuplot)
        )
        .method("draw_domain_pov",
            static_cast<void (container::*)(const char*)>(&container::draw_domain_pov)
        )
        .method(
            "total_particles",
            &container::total_particles
        )
        // Type mismatch from Voro++
        //.method("contains_neighbor", static_cast<void (container::*)(const char*)>(&container::contains_neighbor))
        .method("add_wall!", static_cast<void (container::*)(wall&)>(&container::add_wall))
        .method("add_wall!", static_cast<void (container::*)(wall_list&)>(&container::add_wall))
        .method("point_inside_walls", static_cast<bool (container::*)(double, double, double)>(&container::point_inside_walls))
        // When exporting to Julia, it is needed to defien first voronoicell type
        //.method("apply_walls", static_cast<bool (container::*)(voro::voronoicell&, double, double, double)>(&container::apply_walls))
        ;


    // Class Container Poly
    mod.add_type<container_poly>(
        "RawContainerPoly", jlcxx::julia_type("AbstractRawContainer", "VoroPlusPlus")
    )
        .constructor<double, double, double, double, double, double, int, int, int, bool, bool, bool, int>()
        .method(
            "__cxxwrap_put!",
            static_cast<void (container_poly::*)(int, double, double, double, double)>(&container_poly::put)
        )
        .method("import!", static_cast<void (container_poly::*)(const char*)>(&container_poly::import))
        .method("draw_particles", static_cast<void (container_poly::*)(const char*)>(&container_poly::draw_particles))
        .method("draw_particles_pov", static_cast<void (container_poly::*)(const char*)>(&container_poly::draw_particles_pov))
        .method("draw_cells_gnuplot", static_cast<void (container_poly::*)(const char*)>(&container_poly::draw_cells_gnuplot))
        .method("print_custom", static_cast<void (container_poly::*)(const char*, const char*)>(&container_poly::print_custom))
        .method("draw_cells_pov", static_cast<void (container_poly::*)(const char*)>(&container_poly::draw_cells_pov))
        ;

    // Class Container Iterator
    mod.add_type<c_loop_all>("RawContainerIterator")
        .constructor<container&>()
        .constructor<container_poly&>()
        .method(
            "__cxxwrap_start!",
            &c_loop_all::start
        )
        .method(
            "__cxxwrap_inc!",
            &c_loop_all::inc
        )
        // Inherited from c_loop_base
        .method("ci_x", &c_loop_all::x)
        .method("ci_y", &c_loop_all::y)
        .method("ci_z", &c_loop_all::z)
        .method("ci_pid", &c_loop_all::pid)
        ;
    

    // Class Container Iterator Subset
    mod.add_type<c_loop_subset>("ContainerSubsetIterator")
        .constructor<container&>()
        .constructor<container_poly&>()
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
        .method(
            "__cxxwrap_set_bounds_to_intbox!",
            static_cast<void (c_loop_subset::*)(int, int, int, int, int, int)>(&c_loop_subset::setup_intbox)
        )
        // Inherited from c_loop_base
        .method("cis_x", &c_loop_subset::x)
        .method("cis_y", &c_loop_subset::y)
        .method("cis_z", &c_loop_subset::z)
        .method("cis_pid", &c_loop_subset::pid)
    ;


    // Class Container Iterator Order
    mod.add_type<c_loop_order>("InsertionOrderIterator")
        .constructor<container&, particle_order&>()
        .constructor<container_poly&, particle_order&>()
        .method(
            "__cxxwrap_start!",
            &c_loop_order::start
        )
        .method(
            "__cxxwrap_inc!",
            &c_loop_order::inc
        )
        // Inherited from c_loop_base
        .method("cio_x", &c_loop_order::x)
        .method("cio_y", &c_loop_order::y)
        .method("cio_z", &c_loop_order::z)
        .method("cio_pid", &c_loop_order::pid)
        ;


    // Class VoronoiCell
    mod.add_type<voronoicell_neighbor>(
        "VoronoiCell", jlcxx::julia_type("AbstractVoronoiCell", "VoroPlusPlus")
    )
        .constructor<>()
        .constructor<double>()
        .constructor<container&>()
        .constructor<container_poly&>()
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
            "__cxxwrap_plane!",
            static_cast<bool (voronoicell_neighbor::*)(double, double, double, double)>(&voronoicell_neighbor::plane)
        )
        .method(
            "__cxxwrap_plane!",
            static_cast<bool (voronoicell_neighbor::*)(double, double, double)>(&voronoicell_neighbor::plane)
        )
        .method("volume", &voronoicell_neighbor::volume)
        //.method("check_relations", &voronoicell_neighbor::check_relations)
        //.method("check_duplicates", &voronoicell_neighbor::check_duplicates)
        //.method("max_radius_squared", &voronoicell_neighbor::max_radius_squared)
        .method("number_of_edges", &voronoicell_neighbor::number_of_edges)
        //.method("total_edge_distance", &voronoicell_neighbor::total_edge_distance)
        .method("number_of_faces", &voronoicell_neighbor::number_of_faces)
        .method("surface_area", &voronoicell_neighbor::surface_area)
        .method("draw_gnuplot!", static_cast<void (voronoicell_neighbor::*)(double, double, double, const char*)>(&voronoicell_neighbor::draw_gnuplot))
        .method("draw_pov", static_cast<void (voronoicell_neighbor::*)(double, double, double, const char*)>(&voronoicell_neighbor::draw_pov))
        .method("draw_pov_mesh", static_cast<void (voronoicell_neighbor::*)(double, double, double, const char*)>(&voronoicell_neighbor::draw_pov_mesh))
        .method("plane_intersects", static_cast<bool (voronoicell_neighbor::*)(double, double, double, double)>(&voronoicell_neighbor::plane_intersects))
        .method("nplane", static_cast<bool (voronoicell_neighbor::*)(double, double, double, double, int)>(&voronoicell_neighbor::nplane))
        .method("nplane", static_cast<bool (voronoicell_neighbor::*)(double, double, double, int)>(&voronoicell_neighbor::nplane))
        // Inherited from voronoicell_base
        .method(
            "__cxxwrap_translate!",
            &voronoicell_neighbor::translate
        )
        .method("plane_intersects_guess", static_cast<bool (voronoicell_neighbor::*)(double, double, double, double)>(&voronoicell_neighbor::plane_intersects_guess))
        .method("construct_relations", &voronoicell_neighbor::construct_relations)
        .method("print_edges", &voronoicell_neighbor::print_edges)
        .method("__cycle_up", static_cast<int (voronoicell_neighbor::*)(int, int)>(&voronoicell_neighbor::cycle_up))
        .method("__cycle_down", static_cast<int (voronoicell_neighbor::*)(int, int)>(&voronoicell_neighbor::cycle_down))

        //.method("output_vertices!", static_cast<void (voronoicell::*)(FILE*)>(&voronoicell::output_vertices))
        //.method("output_vertices!", static_cast<void (voronoicell::*)(double, double, double, FILE*)>(&voronoicell::output_vertices));
        //.method("draw_gnuplot", static_cast<void (voronoicell::*)(double, double, double, FILE*)>(&voronoicell::draw_gnuplot))
        //.method("reset_edges", &voronoicell::reset_edges)
        //.method("reset_edges", [](const voronoicell& c) { return dynamic_cast<const C*>(&c)->data; });
        ;

    // Class Containter Periodic Poly (conprdply)
    mod.add_type<container_periodic_poly>("RawContainerTriclinic")
        .constructor<double, double, double, double, double, double, int, int, int, int>()
        .method("__cxxwrap_put!", static_cast<void (container_periodic_poly::*)(int, double, double, double, double)>(&container_periodic_poly::put))
        // Type mismatch from Voro++
        //.method("conprdply_compute_ghost_cell", static_cast<bool (container_periodic_poly::*)(voronoicell&, double, double, double, double)>(&container_periodic_poly::compute_ghost_cell))
        .method("conprdply_draw_particles", static_cast<void (container_periodic_poly::*)(const char*)>(&container_periodic_poly::draw_particles))
        .method("conprdply_draw_cells_gnuplot", static_cast<void (container_periodic_poly::*)(const char*)>(&container_periodic_poly::draw_cells_gnuplot))
        .method("conprdply_draw_domain_gnuplot", static_cast<void (container_periodic_poly::*)(const char*)>(&container_periodic_poly::draw_domain_gnuplot))
        ;
        

    // Class Wall_Sphere
    mod.add_type<wall_sphere>("Wall_Sphere")
        .constructor<double, double, double, double, int>()
        .method("point_inside_wph", static_cast<bool (wall_sphere::*)(double, double, double)>(&wall_sphere::point_inside))
        .method("__cxxwrap_cut_cell!", static_cast<bool (wall_sphere::*)(voronoicell_neighbor&,double,double,double)>(&wall_sphere::cut_cell))
        ;

    // lambdas for container
    mod.method(
        "__cxxwrap_find_cell",
        static_cast<std::tuple<bool,int,double,double,double> (*)(container&, double, double, double)>(
            [] (container &con, double x, double y, double z)
            {
                int pid;
                double rx, ry, rz;
                bool found = con.find_voronoi_cell(x, y, z, rx, ry, rz, pid);
                return std::make_tuple(found, pid, rx, ry, rz);
            }
        )
    );
    mod.method(
        "__cxxwrap_periodic",
        static_cast<std::tuple<bool,bool,bool> (*)(container&)>(
            [] (container &con)
            {
                return std::make_tuple(con.xperiodic, con.yperiodic, con.zperiodic);
            }
        )
    );
    mod.method(
        "__cxxwrap_bounds",
        static_cast<std::tuple<double,double,double,double,double,double> (*)(container&)>(
            [] (container &con)
            {
                return std::make_tuple(con.ax, con.ay, con.az, con.bx, con.by, con.bz);
            }
        )
    );

    // lambdas for voronoicell
    mod.method(
        "centroid",
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
    mod.method(
        "__cxxwrap_get_neighbors!",
        static_cast<void (*)(std::vector<int>&, voronoicell_neighbor&)>(
            [] (std::vector<int>& v, voronoicell_neighbor& vc) { vc.neighbors(v); }   
        )
    );

    // lambdas for loops

    mod.method(
        "pos",
        static_cast<std::tuple<double,double,double> (*)(c_loop_all&)>(
            [] (c_loop_all &cla)
            {
                double x, y, z;
                cla.pos(x, y, z);
                return std::make_tuple(x, y, z);
            }
        )
    );
    mod.method(
        "pos",
        static_cast<std::tuple<double,double,double> (*)(c_loop_order&)>(
            [] (c_loop_order &cla)
            {
                double x, y, z;
                cla.pos(x, y, z);
                return std::make_tuple(x, y, z);
            }
        )
    );
    mod.method(
        "pos",
        static_cast<std::tuple<double,double,double> (*)(c_loop_subset&)>(
            [] (c_loop_subset &cla)
            {
                double x, y, z;
                cla.pos(x, y, z);
                return std::make_tuple(x, y, z);
            }
        )
    );
    mod.method(
        "particle_info",
        static_cast<std::tuple<int,double,double,double,double> (*)(c_loop_all&)>(
            [] (c_loop_all &cla)
            {
                int pid;
                double x, y, z, r;
                cla.pos(pid, x, y, z, r);
                return std::make_tuple(pid, x, y, z, r);
            }
        )
    );
    mod.method(
        "particle_info",
        static_cast<std::tuple<int,double,double,double,double> (*)(c_loop_order&)>(
            [] (c_loop_order &cla)
            {
                int pid;
                double x, y, z, r;
                cla.pos(pid, x, y, z, r);
                return std::make_tuple(pid, x, y, z, r);
            }
        )
    );
    mod.method(
        "particle_info",
        static_cast<std::tuple<int,double,double,double,double> (*)(c_loop_subset&)>(
            [] (c_loop_subset &cla)
            {
                int pid;
                double x, y, z, r;
                cla.pos(pid, x, y, z, r);
                return std::make_tuple(pid, x, y, z, r);
            }
        )
    );

    // Anonymus functions for special cases when two types are needed
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container&, c_loop_all&)>(
            [] (voronoicell_neighbor& vc, container& con, c_loop_all& itr) {
                return con.compute_cell(vc, itr);
            }
        )
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container&, c_loop_order&)>(
            [] (voronoicell_neighbor& vc, container& con, c_loop_order& itr) {
                return con.compute_cell(vc, itr);
            }
        )
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container_poly&, c_loop_all&)>(
            [] (voronoicell_neighbor& vc, container_poly& con, c_loop_all& itr) {
                return con.compute_cell(vc, itr);
            }
        )
    );
    mod.method(
        "__cxxwrap_compute_cell!",
        static_cast<bool (*)(voronoicell_neighbor&, container_poly&, c_loop_order&)>(
            [] (voronoicell_neighbor& vc, container_poly& con, c_loop_order& itr) {
                return con.compute_cell(vc, itr);
            }
        )
    );
    mod.method("__cxxwrap_compute_ghost_cell!", [] (voronoicell_neighbor& vc, container& con, double x, double y, double z) {return con.compute_ghost_cell(vc, x, y, z);});
    mod.method("apply_walls!", [] (voronoicell_neighbor& vc, container& con, double x, double y, double z){ return con.apply_walls(vc, x, y, z);});

    // Public Menbers from Container Class
    mod.method("get_particle_id", &get_particle_id);

    ///////////////////////// refactors for Ref substitution //////////////////////

    // mod.method("find_voro_cell", [] (container& con, double x, double y, double z) { double rx, ry, rz; 
    //     int pid; bool found = con.find_voronoi_cell(x, y, z, rx, ry, rz, pid); 
    //     return std::make_tuple(found, rx, ry, rz, pid);});

    // mod.method("get_centroid", [] (voronoicell& v) {double x, y, z; v.centroid(x, y, z); return std::make_tuple(x,y,z);});

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}