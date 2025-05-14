#include "jlcxx/jlcxx.hpp"
#include "jlcxx/tuple.hpp"
#include "voro++.hh"


/* Get and Set methods for accessing public members from voronoicell class cel.hh*/

int get_up(voro::voronoicell& v) {

    return v.up;
}

int get_p(voro::voronoicell& v) {

    return v.p;
}

int get_current_vertices(voro::voronoicell& v) {

    return v.current_vertices;
}

int* get_nu(voro::voronoicell& v) {

    return v.nu;
}

int get_edge(voro::voronoicell& v, int i, int j) {

    return v.ed[i][j];
}

int set_edge(voro::voronoicell& v, int k, int i, int j) {

    v.ed[i][j] = k;
    return k;
}

double* get_pts(voro::voronoicell& v) {

    return v.pts;
}


/* Get and Set methods for accessing public members from container class container.hh */
int get_particle_id(voro::container& con, int i, int j) {

    return con.id[i][j];
}



/* Extern Method for special cases Voronoicell */
extern "C" {

    void draw_gnuplot_voronoicell(voro::voronoicell* vc, double x, double y, double z, FILE *fp) { vc->draw_gnuplot(x, y, z, fp); }

    // Version implemented from Voro++  >>>  void output_vertices(FILE *fp=stdout), ignoring stdout by default
    // without position parameters
    void output_vertices_nopos(voro::voronoicell* vc, FILE *fp) { vc->output_vertices(fp); }

    // Version implemented from Voro++  >>>  void output_vertices(double x,double y,double z,FILE *fp=stdout), ignoring stdout by default
    // with position parameters
    void output_vertices_positions(voro::voronoicell* vc, double x, double y, double z, FILE *fp) { vc->output_vertices(x, y, z, fp); }

    // Version implemented from Voro++  >>>  void output_vertex_orders(FILE *fp=stdout), ignoring stdout by default
    void output_vertex_orders_vorocell(voro::voronoicell* vc, FILE *fp) { vc->output_vertex_orders(fp); }

    // Version implemented from Voro++  >>>  inline void output_face_perimeters(FILE *fp=stdout), ignoring stdout by default
    void output_face_perimeters_vorocell(voro::voronoicell* vc, FILE *fp) { vc->output_face_perimeters(fp); }

    // Version implemented from Voro++  >>>  inline void output_face_freq_table(FILE *fp=stdout), ignoring stdout by default
    void output_face_freq_table_vorocell(voro::voronoicell* vc, FILE *fp) { vc->output_face_freq_table(fp); }

    // Version implemented from Voro++  >>>  inline void output_face_orders(FILE *fp=stdout), ignoring stdout by default
    void output_face_orders_vorocell(voro::voronoicell* vc, FILE *fp) { vc->output_face_orders(fp); }
    
    // Version implemented from Voro++  >>>  inline void output_face_areas(FILE *fp=stdout), ignoring stdout by default
    void output_face_areas_vorocell(voro::voronoicell* vc, FILE *fp) { vc->output_face_areas(fp); }
    
    // Version implemented from Voro++  >>>  inline void output_normals(FILE *fp=stdout), ignoring stdout by default
    void output_normals_vorocell(voro::voronoicell* vc, FILE *fp) { vc->output_normals(fp); }
    
    // Version implemented from Voro++  >>>  inline void output_face_vertices(FILE *fp=stdout), ignoring stdout by default
    void output_face_vertices_vorocell(voro::voronoicell* vc, FILE *fp) { vc->output_face_vertices(fp); }

    bool compute_ghost_cell_conprdply(voro::container_periodic_poly* con, voro::voronoicell* c, double x, double y, double z, double r) {

        return con->compute_ghost_cell(c, x, y, z, r);
    }

}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    using namespace voro;

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
        "RawContainer", jlcxx::julia_type("AbstractContainer", "VoroPlusPlus")
    )
        .constructor<double, double, double, double, double, double, int, int, int, bool, bool, bool, int>()
        .method(
            "__cxxwrap_add_point!",
            static_cast<void (container::*)(int, double, double, double)>(&container::put)
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
        .method("__cxxwrap_find_cell",
            [] (container &con, double x, double y, double z)
            {
                int pid;
                double rx, ry, rz;
                bool found = con.find_voronoi_cell(x, y, z, rx, ry, rz, pid);
                return std::make_tuple(found, pid, rx, ry, rz);
            }
        )
        .method(
            "clear!",
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
        .method("add_wall!", static_cast<void (container::*)(wall*)>(&container::add_wall))
        .method("add_wall!", static_cast<void (container::*)(wall&)>(&container::add_wall))
        .method("add_wall!", static_cast<void (container::*)(wall_list&)>(&container::add_wall))
        .method("point_inside_walls", static_cast<bool (container::*)(double, double, double)>(&container::point_inside_walls))
        // When exporting to Julia, it is needed to defien first voronoicell type
        //.method("apply_walls", static_cast<bool (container::*)(voro::voronoicell&, double, double, double)>(&container::apply_walls))
        ;


    // Class Container Poly
    mod.add_type<container_poly>(
        "RawContainerPoly", jlcxx::julia_type("AbstractContainer", "VoroPlusPlus")
    )
        .constructor<double, double, double, double, double, double, int, int, int, bool, bool, bool, int>()
        .method(
            "__cxxwrap_add_point!",
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
        .method("__cxxwrap_start!", &c_loop_all::start)
        .method("__cxxwrap_next!", &c_loop_all::inc)
        // Inherited from c_loop_base
        .method("pos", [] (c_loop_all &cla, double &x, double &y, double &z){cla.pos(x, y, z);})
        .method("pos", [] (c_loop_all &cla, int &pid, double &x, double &y, double &z, double &r){cla.pos(pid, x, y, z, r);})
        .method("ci_x", &c_loop_all::x)
        .method("ci_y", &c_loop_all::y)
        .method("ci_z", &c_loop_all::z)
        .method("ci_pid", &c_loop_all::pid)
        ;
    

    // Class Container Iterator Subset
    mod.add_type<c_loop_subset>("Container_Iterator_Subset")
        .constructor<container&>()
        .constructor<container_poly&>()
        .method("cis_start", &c_loop_subset::start)
        .method("cis_next", &c_loop_subset::inc)
        .method("cis_setup_sphere", static_cast<void (c_loop_subset::*)(double, double, double, double, bool)>(&c_loop_subset::setup_sphere))
        .method("cis_setup_box", static_cast<void (c_loop_subset::*)(double, double, double, double, double, double, bool)>(&c_loop_subset::setup_box))
        .method("cis_setup_intbox", static_cast<void (c_loop_subset::*)(int, int, int, int, int, int)>(&c_loop_subset::setup_intbox))
        // Inherited from c_loop_base
        .method("cis_pos", [] (c_loop_subset &cls, double &x, double &y, double &z){cls.pos(x, y, z);})
        .method("cis_pos", [] (c_loop_subset &cls, int &pid, double &x, double &y, double &z, double &r){cls.pos(pid, x, y, z, r);})
        .method("cis_x", &c_loop_subset::x)
        .method("cis_y", &c_loop_subset::y)
        .method("cis_z", &c_loop_subset::z)
        .method("cis_pid", &c_loop_subset::pid)
        ;
    
    
    // Class Particle Order
    mod.add_type<particle_order>("Particle_Order")
        .constructor<int>()
        .method("po_add", static_cast<void (particle_order::*)(int, int)>(&particle_order::add))
        ;


    // Class Container Iterator Order
    mod.add_type<c_loop_order>("Container_Iterator_Order")
        .constructor<container&, particle_order&>()
        .constructor<container_poly&, particle_order&>()
        .method("cio_start", &c_loop_order::start)
        .method("cio_next", &c_loop_order::inc)
        // Inherited from c_loop_base
        .method("cio_pos", [] (c_loop_order &clo, double &x, double &y, double &z){clo.pos(x, y, z);})
        .method("cio_pos", [] (c_loop_order &clo, int &pid, double &x, double &y, double &z, double &r){clo.pos(pid, x, y, z, r);})
        .method("cio_x", &c_loop_order::x)
        .method("cio_y", &c_loop_order::y)
        .method("cio_z", &c_loop_order::z)
        .method("cio_pid", &c_loop_order::pid)
        ;


    // Class VoronoiCell
    mod.add_type<voronoicell>("VoronoiCell")
        .constructor<>()
        .constructor<double>()
        .constructor<container&>()
        .method("init!", &voronoicell::init)
        .method("init_l_shape!", &voronoicell::init_l_shape)
        .method("add_plane!", static_cast<bool (voronoicell::*)(double, double, double, double)>(&voronoicell::plane))
        .method("add_plane!", static_cast<bool (voronoicell::*)(double, double, double)>(&voronoicell::plane))
        .method("volume", &voronoicell::volume)
        .method("check_relations", &voronoicell::check_relations)
        .method("check_duplicates", &voronoicell::check_duplicates)
        .method("max_radius_squared", &voronoicell::max_radius_squared)
        .method("number_of_edges", &voronoicell::number_of_edges)
        .method("total_edge_distance", &voronoicell::total_edge_distance)
        .method("number_of_faces", &voronoicell::number_of_faces)
        .method("surface_area", &voronoicell::surface_area)
        .method("draw_gnuplot!", static_cast<void (voronoicell::*)(double, double, double, const char*)>(&voronoicell::draw_gnuplot))
        .method("draw_pov", static_cast<void (voronoicell::*)(double, double, double, const char*)>(&voronoicell::draw_pov))
        .method("draw_pov_mesh", static_cast<void (voronoicell::*)(double, double, double, const char*)>(&voronoicell::draw_pov_mesh))
        .method("init_octahedron", static_cast<void (voronoicell::*)(double)>(&voronoicell::init_octahedron))
        .method("plane_intersects", static_cast<bool (voronoicell::*)(double, double, double, double)>(&voronoicell::plane_intersects))
        .method("centroid!", static_cast<void (voronoicell::*)(double&, double&, double&)>(&voronoicell::centroid))

        .method("nplane", static_cast<bool (voronoicell::*)(double, double, double, double, int)>(&voronoicell::nplane))
        .method("nplane", static_cast<bool (voronoicell::*)(double, double, double, int)>(&voronoicell::nplane))
        .method("init_tetrahedron", static_cast<void (voronoicell::*)(double, double, double, double, double, double, double, double, double, double, double, double)>(&voronoicell::init_tetrahedron))

        // Inherited from voronoicell_base
        .method("init_base", static_cast<void (voronoicell::*)(double, double, double, double, double, double)>(&voronoicell::init_base))
        .method("init_octahedron_base", static_cast<void (voronoicell::*)(double)>(&voronoicell::init_octahedron_base))
        .method("init_tetrahedron_base", static_cast<void (voronoicell::*)(double, double, double, double, double, double, double, double, double, double, double, double)>(&voronoicell::init_tetrahedron_base))
        .method("translate", static_cast<void (voronoicell::*)(double, double, double)>(&voronoicell::translate))
        .method("plane_intersects_guess", static_cast<bool (voronoicell::*)(double, double, double, double)>(&voronoicell::plane_intersects_guess))
        .method("construct_relations", &voronoicell::construct_relations)
        .method("print_edges", &voronoicell::print_edges)
        .method("cycle_up", static_cast<int (voronoicell::*)(int, int)>(&voronoicell::cycle_up))
        .method("cycle_down", static_cast<int (voronoicell::*)(int, int)>(&voronoicell::cycle_down))

        //.method("output_vertices!", static_cast<void (voronoicell::*)(FILE*)>(&voronoicell::output_vertices))
        //.method("output_vertices!", static_cast<void (voronoicell::*)(double, double, double, FILE*)>(&voronoicell::output_vertices));
        //.method("draw_gnuplot", static_cast<void (voronoicell::*)(double, double, double, FILE*)>(&voronoicell::draw_gnuplot))
        //.method("reset_edges", &voronoicell::reset_edges)
        //.method("reset_edges", [](const voronoicell& c) { return dynamic_cast<const C*>(&c)->data; });
        ;

    
    // Class VoronoiCell_Neighbor
    mod.add_type<voronoicell_neighbor>("VoronoiCell_Neighbor")
        .constructor<>()
        .constructor<double>()
        .constructor<container&>()
        .method("init_vcn", static_cast<void (voronoicell_neighbor::*)(double,double,double,double,double,double)>(&voronoicell_neighbor::init))
        .method("init_octahedron_vcn", static_cast<void (voronoicell_neighbor::*)(double)>(&voronoicell_neighbor::init_octahedron))
        .method("init_tetrahedron_vcn", static_cast<void (voronoicell_neighbor::*)(double,double,double,double,double,double,double,double,double,double,double,double)>(&voronoicell_neighbor::init_tetrahedron))
        .method("nplane_rsq_vcn", static_cast<bool (voronoicell_neighbor::*)(double,double,double,double,int)>(&voronoicell_neighbor::nplane))
        .method("nplane_vcn", static_cast<bool (voronoicell_neighbor::*)(double,double,double,int)>(&voronoicell_neighbor::nplane))
        .method("plane_rsq_vcn", static_cast<bool (voronoicell_neighbor::*)(double,double,double,double)>(&voronoicell_neighbor::plane))
        .method("plane_vcn", static_cast<bool (voronoicell_neighbor::*)(double,double,double)>(&voronoicell_neighbor::plane))
        .method("check_facets_vcn", &voronoicell_neighbor::check_facets)
        .method("print_edges_neighbors_vcn", static_cast<void (voronoicell_neighbor::*)(int)>(&voronoicell_neighbor::print_edges_neighbors))
        .method("neighbors_vcn", static_cast<void (voronoicell_neighbor::*)(std::vector<int>&)>(&voronoicell_neighbor::neighbors))
        // Inherited from voronoicel_base
        .method("translate_vcn", static_cast<void (voronoicell_neighbor::*)(double, double, double)>(&voronoicell_neighbor::translate))
        .method("draw_pov_vcn", static_cast<void (voronoicell_neighbor::*)(double, double, double, const char*)>(&voronoicell_neighbor::draw_pov))
        .method("draw_pov_mesh_vcn", static_cast<void (voronoicell_neighbor::*)(double, double, double, const char*)>(&voronoicell_neighbor::draw_pov_mesh))
        .method("draw_gnuplot_vcn", static_cast<void (voronoicell_neighbor::*)(double, double, double, const char*)>(&voronoicell_neighbor::draw_gnuplot))
        .method("volume_vcn", &voronoicell_neighbor::volume)
        .method("max_radius_squared_vcn", &voronoicell_neighbor::max_radius_squared)
        .method("total_edge_distance_vcn", &voronoicell_neighbor::total_edge_distance)
        .method("surface_area_vcn", &voronoicell_neighbor::surface_area)
        .method("centroid_vcn", static_cast<void (voronoicell_neighbor::*)(double&, double&, double&)>(&voronoicell_neighbor::centroid))
        .method("number_of_faces_vcn", &voronoicell_neighbor::number_of_faces)
        .method("number_of_edges_vcn", &voronoicell_neighbor::number_of_edges)
        ;


    // Class Containter Periodic Poly (conprdply)
    mod.add_type<container_periodic_poly>("Container_Periodic_Poly")
        .constructor<double, double, double, double, double, double, int, int, int, int>()
        .method("conprdply_add_point!", static_cast<void (container_periodic_poly::*)(int, double, double, double, double)>(&container_periodic_poly::put))
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
        .method("cut_cell_vc", static_cast<bool (wall_sphere::*)(voronoicell&,double,double,double)>(&wall_sphere::cut_cell))
        .method("cut_cell_vcn", static_cast<bool (wall_sphere::*)(voronoicell_neighbor&,double,double,double)>(&wall_sphere::cut_cell))
        ;

    
    // Public Members from VoronoiCell Class
    mod.method("root_vertex", &get_up);
    mod.method("num_vertices", &get_p);
    mod.method("__get_nu", &get_nu);
    mod.method("__get_current_vertices", &get_current_vertices);
    mod.method("__get_edge", &get_edge);
    mod.method("__set_edge!", &set_edge);
    mod.method("__get_pts", &get_pts);

    // Anonymus functions for special cases when two types are needed
    mod.method("compute_cell!", [] (voronoicell& vc, container& con, c_loop_all& itr) {return con.compute_cell(vc, itr);});
    mod.method("compute_cell!", [] (voronoicell& vc, container& con, int ijk, int q) {return con.compute_cell(vc, ijk, q);});
    mod.method("compute_ghost_cell!", [] (voronoicell& vc, container& con, double x, double y, double z) {return con.compute_ghost_cell(vc, x, y, z);});
    mod.method("apply_walls!", [] (voronoicell& vc, container& con, double x, double y, double z){ return con.apply_walls(vc, x, y, z);});

    // Public Menbers from Container Class
    mod.method("get_particle_id", &get_particle_id);

    ///////////////////////// refactors for Ref substitution //////////////////////

    mod.method("get_pos", [] (c_loop_all& cla) {int pid; double x, y, z, r; cla.pos(pid, x, y, z, r); return std::make_tuple(pid,x,y,z,r);});

    mod.method("find_voro_cell", [] (container& con, double x, double y, double z) { double rx, ry, rz; 
        int pid; bool found = con.find_voronoi_cell(x, y, z, rx, ry, rz, pid); 
        return std::make_tuple(found, rx, ry, rz, pid);});

    mod.method("get_centroid", [] (voronoicell& v) {double x, y, z; v.centroid(x, y, z); return std::make_tuple(x,y,z);});

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}