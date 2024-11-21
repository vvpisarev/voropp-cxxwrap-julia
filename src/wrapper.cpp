#include "jlcxx/jlcxx.hpp"
#include "voro++.hh"

int get_up(voro::voronoicell& v)
{
    return v.up;
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    using namespace voro;

    mod.add_type<container>("Container")
        .constructor<double, double, double, double, double, double, int, int, int, bool, bool, bool, int>()
        .method("add_point!", static_cast<void (container::*)(int, double, double, double)>(&container::put))
        .method("import!", static_cast<void (container::*)(const char*)>(&container::import))
        //.method("draw_cells_gnuplot", static_cast<void (container::*)(FILE*)>(&container::draw_cells_gnuplot))
        .method("draw_particles", static_cast<void (container::*)(const char*)>(&container::draw_particles))
        .method("draw_cells_gnuplot", static_cast<void (container::*)(const char*)>(&container::draw_cells_gnuplot));

    mod.add_type<c_loop_all>("ContainerIterator")
        .constructor<container&>()
        .method("start!", &c_loop_all::start)
        .method("next!", &c_loop_all::inc);

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
        //.method("draw_gnuplot", static_cast<void (voronoicell::*)(double, double, double, FILE*)>(&voronoicell::draw_gnuplot))
        .method("draw_gnuplot", static_cast<void (voronoicell::*)(double, double, double, const char*)>(&voronoicell::draw_gnuplot));
    mod.method("get_up", &get_up);

    mod.method("compute_cell!", [] (voronoicell& vc, container& con, c_loop_all& itr) {return con.compute_cell(vc, itr);});
    mod.method("compute_cell!", [] (voronoicell& vc, container& con, int ijk, int q) {return con.compute_cell(vc, ijk, q);});
    
}
