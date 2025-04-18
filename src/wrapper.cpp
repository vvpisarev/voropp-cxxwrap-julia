#include "jlcxx/jlcxx.hpp"
#include "voro++.hh"


/* Get and Set methods for accessing public members from voronoicell class cel.hh*/

int get_up(voro::voronoicell& v)
{
    return v.up;
}

int get_p(voro::voronoicell& v)
{
    return v.p;
}

int get_current_vertices(voro::voronoicell& v)
{
    return v.current_vertices;
}

int* get_nu(voro::voronoicell& v)
{
    return v.nu;
}

int get_edge(voro::voronoicell& v, int i, int j)
{
    return v.ed[i][j];
}

int set_edge(voro::voronoicell& v, int k, int i, int j)
{
    v.ed[i][j] = k;
    return k;
}

double* get_pts(voro::voronoicell& v)
{
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
    mod.add_type<container>("Container")
        .constructor<double, double, double, double, double, double, int, int, int, bool, bool, bool, int>()
        .method("add_point!", static_cast<void (container::*)(int, double, double, double)>(&container::put))
        .method("import!", static_cast<void (container::*)(const char*)>(&container::import))
        // When exporting to Julia, it is no possible to derive automatically type FILE
        //.method("draw_cells_gnuplot", static_cast<void (container::*)(FILE*)>(&container::draw_cells_gnuplot))
        .method("draw_particles", static_cast<void (container::*)(const char*)>(&container::draw_particles))
        .method("draw_particles_pov", static_cast<void (container::*)(const char*)>(&container::draw_particles_pov))
        .method("draw_cells_gnuplot", static_cast<void (container::*)(const char*)>(&container::draw_cells_gnuplot))
        .method("print_custom!", static_cast<void (container::*)(const char*, const char*)>(&container::print_custom))
        .method("draw_cells_pov!", static_cast<void (container::*)(const char*)>(&container::draw_cells_pov))
        .method("find_voronoi_cell", [] (container &con, double x, double y, double z, double &rx, double &ry, double &rz, int &pid){return con.find_voronoi_cell(x, y, z, rx, ry, rz, pid);})
        .method("clear", &container::clear)
        .method("compute_all_cells", &container::compute_all_cells)
        .method("sum_cell_volumes", &container::sum_cell_volumes)
        .method("point_inside!", static_cast<bool (container::*)(double, double, double)>(&container::point_inside))
        .method("region_count", &container::region_count)
        // When exporting to Julia, it is needed to defien first voronoicell type
        //.method("initialize_voronoicell", static_cast<bool (container::*)(voro::voronoicell&, int, int, int, int, int, int&, int&, int&, double&, double&, double&, int&)>(&container::initialize_voronoicell))
        .method("initialize_search", static_cast<void (container::*)(int, int, int, int, int&, int&, int&, int&)>(&container::initialize_search))
        .method("frac_pos", static_cast<void (container::*)(double, double, double, double, double, double, double&, double&, double&)>(&container::frac_pos))
        .method("region_index", static_cast<int (container::*)(int, int, int, int, int, int, double&, double&, double&, int&)>(&container::region_index))
        .method("draw_domain_gnuplot", static_cast<void (container::*)(const char*)>(&container::draw_domain_gnuplot))
        .method("draw_domain_pov", static_cast<void (container::*)(const char*)>(&container::draw_domain_pov))
        .method("total_particles", &container::total_particles)
        // Type mismatch from Voro++
        //.method("contains_neighbor", static_cast<void (container::*)(const char*)>(&container::contains_neighbor))
        .method("add_wall!", static_cast<void (container::*)(wall*)>(&container::add_wall))
        .method("add_wall!", static_cast<void (container::*)(wall&)>(&container::add_wall))
        .method("add_wall!", static_cast<void (container::*)(wall_list&)>(&container::add_wall))
        .method("point_inside_walls", static_cast<bool (container::*)(double, double, double)>(&container::point_inside_walls))
        // When exporting to Julia, it is needed to defien first voronoicell type
        //.method("apply_walls", static_cast<bool (container::*)(voro::voronoicell&, double, double, double)>(&container::apply_walls))
        
        // Protected Methods cannot be added
        //.method("add_particle_memory", static_cast<void (container::*)(int)>(&container::add_particle_memory))
        //.method("put_locate_block", static_cast<bool (container::*)(int&, double&, double&, double&)>(&container::put_locate_block))
        //.method("put_remap", static_cast<bool (container::*)(int&, double&, double&, double&)>(&container::put_remap))
        //.method("remap", static_cast<bool (container::*)(int&, int&, int&, int&, int&, int&, double&, double&, double&, int&)>(&container::remap))
        ;


    // Class Container Poly
    mod.add_type<container_poly>("Container_Poly")
        .constructor<double, double, double, double, double, double, int, int, int, bool, bool, bool, int>()
        .method("conp_add_point!", static_cast<void (container_poly::*)(int, double, double, double, double)>(&container_poly::put))
        .method("conp_import!", static_cast<void (container_poly::*)(const char*)>(&container_poly::import))
        .method("conp_draw_particles", static_cast<void (container_poly::*)(const char*)>(&container_poly::draw_particles))
        .method("conp_draw_particles_pov", static_cast<void (container_poly::*)(const char*)>(&container_poly::draw_particles_pov))
        .method("conp_draw_cells_gnuplot", static_cast<void (container_poly::*)(const char*)>(&container_poly::draw_cells_gnuplot))
        .method("conp_print_custom!", static_cast<void (container_poly::*)(const char*, const char*)>(&container_poly::print_custom))
        .method("conp_draw_cells_pov!", static_cast<void (container_poly::*)(const char*)>(&container_poly::draw_cells_pov));


    mod.add_type<c_loop_all>("Container_Iterator")
        .constructor<container&>()
        .constructor<container_poly&>()
        .method("start!", &c_loop_all::start)
        .method("next!", &c_loop_all::inc)
        .method("pos", [] (c_loop_all &cla, int &pid, double &x, double &y, double &z, double &r){cla.pos(pid, x, y, z, r);});


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

    mod.method("compute_cell!", [] (voronoicell& vc, container& con, c_loop_all& itr) {return con.compute_cell(vc, itr);});
    mod.method("compute_cell!", [] (voronoicell& vc, container& con, int ijk, int q) {return con.compute_cell(vc, ijk, q);});
    mod.method("compute_ghost_cell!", [] (voronoicell& vc, container& con, double x, double y, double z) {return con.compute_ghost_cell(vc, x, y, z);});

    // Public Menbers from Container Class

    mod.method("get_particle_id", &get_particle_id);
}