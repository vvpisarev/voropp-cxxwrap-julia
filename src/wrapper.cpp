#include "jlcxx/jlcxx.hpp"
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

struct box_bounds
{
    double ax;
    double ay;
    double az;
    double bx;
    double by;
    double bz;
};

struct pbc
{
    bool px;
    bool py;
    bool pz;
};

struct vec3d
{
    double x;
    double y;
    double z;
};

struct particle_info
{
    double x;
    double y;
    double z;
    double r;
    int pid;
};

struct loop_indices
{
    int i;
    int j;
    int k;
    int ijk;
    int q;
};

#include "containers.cpp"
#include "loops.cpp"
#include "cell.cpp"

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    using namespace voro;

    mod.map_type<box_bounds>("CxxBoxBounds");
    mod.map_type<pbc>("CxxPBC");
    mod.map_type<vec3d>("CxxVec3D");
    mod.map_type<particle_info>("CxxParticleInfo");
    mod.map_type<loop_indices>("CxxLoopIndices");
    // Class Container
    auto class_container = mod.add_type<container>(
        "Container", jlcxx::julia_type("AbstractContainer")
    )
        .constructor<double, double, double, double, double, double, int, int, int, bool, bool, bool, int>()
    ;

    // Class Container Poly
    auto class_container_poly = mod.add_type<container_poly>(
        "ContainerPoly", jlcxx::julia_type("AbstractContainer")
    )
        .constructor<double, double, double, double, double, double, int, int, int, bool, bool, bool, int>()
    ;

    // Class Containter Periodic
    auto class_container_periodic = mod.add_type<container_periodic>(
        "ContainerTriclinic", jlcxx::julia_type("AbstractContainer")
    )
        .constructor<double, double, double, double, double, double, int, int, int, int>()
    ;

    auto class_container_periodic_poly = mod.add_type<container_periodic_poly>(
        "ContainerTriclinicPoly", jlcxx::julia_type("AbstractContainer")
    )
        .constructor<double, double, double, double, double, double, int, int, int, int>()
    ;

    // Class VoronoiCell
    auto class_voronoicell_neighbor = mod.add_type<voronoicell_neighbor>(
        "VoronoiCell", jlcxx::julia_type("AbstractVoronoiCell")
    )
        .constructor<>()
        .constructor<double>()
        .constructor<container&>()
        .constructor<container_poly&>()
    ;

    // Class Particle Order
    auto class_particle_order = mod.add_type<particle_order>(
        "InsertionOrder",  jlcxx::julia_type("ContainerIterationOrder")
    )
        .constructor<>()
    ;

    // Class Container Iterator
    auto class_c_loop_all = mod.add_type<c_loop_all>("ContainerIterator")
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

    auto class_c_loop_all_periodic = mod.add_type<c_loop_all_periodic>("TriclinicContainerIterator")
        .constructor<container_periodic&>()
        .constructor<container_periodic_poly&>()
    ;

    auto class_c_loop_order_periodic = mod.add_type<c_loop_order_periodic>("TriclinicInsertionOrderIterator")
        .constructor<container_periodic&, particle_order&>()
        .constructor<container_periodic_poly&, particle_order&>()
    ;

    // Class Wall   
    auto class_wall_sphere = mod.add_type<wall_sphere>(
        "WallSphere", jlcxx::julia_type("AbstractWall")
    )
        .constructor<double, double, double, double>()
        .constructor<double, double, double, double, int>()
    ;

    auto class_wall_cylinder = mod.add_type<wall_cylinder>(
        "WallCylinder", jlcxx::julia_type("AbstractWall")
    )
        .constructor<double, double, double, double, double, double, double, int>()
        .constructor<double, double, double, double, double, double, double>()
    ;

    auto class_wall_cone = mod.add_type<wall_cone>(
        "WallCone", jlcxx::julia_type("AbstractWall")
    )
        .constructor<double, double, double, double, double, double, double, int>()
        .constructor<double, double, double, double, double, double, double>()
    ;

    auto class_wall_plane = mod.add_type<wall_plane>(
        "WallPlane", jlcxx::julia_type("AbstractWall")
    )
        .constructor<double, double, double, double, int>()
        .constructor<double, double, double, double>()
    ;

    export_containers_methods(mod);
    export_loops_methods(mod);
    export_voronoicell_methods(mod);

    // Class container_periodic service methods
    // not performance-critical, exporting as std::function
    class_container_periodic
        .method(
            "__cxxwrap_print_all_particles",
            &container_periodic::print_all_particles
        )
        .method(
            "__cxxwrap_check_compartmentalized",
            &container_periodic::check_compartmentalized
        )
        ;
    // Class VoronoiCell service methods
    // not performance-critical, exporting as std::function
    class_voronoicell_neighbor
        .method(
            "__cxxwrap_plane_intersects",
            &voronoicell_neighbor::plane_intersects
        )
        .method(
            "__cxxwrap_plane_intersects_guess",
            &voronoicell_base::plane_intersects_guess
        )
        .method(
            "__cxxwrap_construct_relations!",
            &voronoicell_neighbor::construct_relations
        )
        .method(
            "__cxxwrap_check_relations",
            &voronoicell_neighbor::check_relations
        )
        .method(
            "__cxxwrap_check_duplicates",
            &voronoicell_neighbor::check_duplicates
        )
        .method(
            "__cxxwrap_print_edges",
            &voronoicell_neighbor::print_edges
        )
        .method(
            "__cxxwrap_print_edges_neighbors",
            &voronoicell_neighbor::print_edges_neighbors
        )
        .method(
            "__cxxwrap_check_facets",
            &voronoicell_neighbor::check_facets
        )
    ;

}