void export_loops_methods(jlcxx::Module& mod)
{
    using namespace voro;

    // lambdas for loop classes

    auto __cxxwrap_start = [] (auto& cl) { return cl.start(); };
    auto __cxxwrap_inc = [] (auto& cl) { return cl.inc(); };
    auto __cxxwrap_pos = [] (auto &cl) 
        {
            double x, y, z;
            cl.pos(x, y, z);
            return std::make_tuple(x, y, z);
        };

    auto __cxxwrap_particle_info = [] (auto &cl)
        {
            int pid;
            double x, y, z, r;
            cl.pos(pid, x, y, z, r);
            return std::make_tuple(pid, x, y, z, r);
        };

    auto __cxxwrap_loop_state = [] (auto &cl)
        {
            return std::make_tuple(cl.i, cl.j, cl.k, cl.ijk, cl.q);
        };

    auto __cxxwrap_restore_loop_state = [] (auto &cl, int i, int j, int k, int ijk, int q)
        {
            cl.i = i;
            cl.j = j;
            cl.k = k;
            cl.ijk = ijk;
            cl.q = q;
        };

    mod.method(
        "__cxxwrap_start!",
        static_cast<bool (*)(c_loop_all&)>(
            __cxxwrap_start
        )
    );
    mod.method(
        "__cxxwrap_inc!",
        static_cast<bool (*)(c_loop_all&)>(
            __cxxwrap_inc
        )
    );
    mod.method(
        "__cxxwrap_pos",
        static_cast<std::tuple<double,double,double> (*)(c_loop_all&)>(
            __cxxwrap_pos
        )
    );
    mod.method(
        "__cxxwrap_particle_info",
        static_cast<std::tuple<int,double,double,double,double> (*)(c_loop_all&)>(
            __cxxwrap_particle_info
        )
    );
    mod.method(
        "__cxxwrap_start!",
        static_cast<bool (*)(c_loop_order&)>(
            __cxxwrap_start
        )
    );
    mod.method(
        "__cxxwrap_inc!",
        static_cast<bool (*)(c_loop_order&)>(
            __cxxwrap_inc
        )
    );
    mod.method(
        "__cxxwrap_pos",
        static_cast<std::tuple<double,double,double> (*)(c_loop_order&)>(
            __cxxwrap_pos
        )
    );
    mod.method(
        "__cxxwrap_particle_info",
        static_cast<std::tuple<int,double,double,double,double> (*)(c_loop_order&)>(
            __cxxwrap_particle_info
        )
    );
    mod.method(
        "__cxxwrap_start!",
        static_cast<bool (*)(c_loop_subset&)>(
            __cxxwrap_start
        )
    );
    mod.method(
        "__cxxwrap_inc!",
        static_cast<bool (*)(c_loop_subset&)>(
            __cxxwrap_inc
        )
    );
    mod.method(
        "__cxxwrap_pos",
        static_cast<std::tuple<double,double,double> (*)(c_loop_subset&)>(
            __cxxwrap_pos
        )
    );
    mod.method(
        "__cxxwrap_particle_info",
        static_cast<std::tuple<int,double,double,double,double> (*)(c_loop_subset&)>(
            __cxxwrap_particle_info
        )
    );
    mod.method(
        "__cxxwrap_setup_sphere!",
        static_cast<void (*)(c_loop_subset&, double, double, double, double, bool)>(
            [](c_loop_subset& cl, double x, double y, double z, double r, bool bounds_test)
            {
                cl.setup_sphere(x, y, z, r, bounds_test);
            }
        )
    );
    mod.method(
        "__cxxwrap_setup_box!",
        static_cast<void (*)(c_loop_subset&, double, double, double, double, double, double, bool)>(
            [](c_loop_subset& cl, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, bool bounds_test)
            {
                cl.setup_box(xmin, xmax, ymin, ymax, zmin, zmax, bounds_test);
            }
        )
    );

    mod.method(
        "__cxxwrap_start!",
        static_cast<bool (*)(c_loop_all_periodic&)>(
            __cxxwrap_start
        )
    );
    mod.method(
        "__cxxwrap_inc!",
        static_cast<bool (*)(c_loop_all_periodic&)>(
            __cxxwrap_inc
        )
    );
    mod.method(
        "__cxxwrap_pos",
        static_cast<std::tuple<double,double,double> (*)(c_loop_all_periodic&)>(
            __cxxwrap_pos
        )
    );
    mod.method(
        "__cxxwrap_particle_info",
        static_cast<std::tuple<int,double,double,double,double> (*)(c_loop_all_periodic&)>(
            __cxxwrap_particle_info
        )
    );
    mod.method(
        "__cxxwrap_start!",
        static_cast<bool (*)(c_loop_order_periodic&)>(
            __cxxwrap_start
        )
    );
    mod.method(
        "__cxxwrap_inc!",
        static_cast<bool (*)(c_loop_order_periodic&)>(
            __cxxwrap_inc
        )
    );
    mod.method(
        "__cxxwrap_pos",
        static_cast<std::tuple<double,double,double> (*)(c_loop_order_periodic&)>(
            __cxxwrap_pos
        )
    );
    mod.method(
        "__cxxwrap_particle_info",
        static_cast<std::tuple<int,double,double,double,double> (*)(c_loop_order_periodic&)>(
            __cxxwrap_particle_info
        )
    );

    mod.method(
        "__cxxwrap_loop_state",
        static_cast<std::tuple<int,int,int,int,int> (*)(c_loop_all&)>(
            __cxxwrap_loop_state
        )
    );
    mod.method(
        "__cxxwrap_loop_state",
        static_cast<std::tuple<int,int,int,int,int> (*)(c_loop_order&)>(
            __cxxwrap_loop_state
        )
    );
    mod.method(
        "__cxxwrap_loop_state",
        static_cast<std::tuple<int,int,int,int,int> (*)(c_loop_subset&)>(
            __cxxwrap_loop_state
        )
    );
    mod.method(
        "__cxxwrap_loop_state",
        static_cast<std::tuple<int,int,int,int,int> (*)(c_loop_all_periodic&)>(
            __cxxwrap_loop_state
        )
    );
    mod.method(
        "__cxxwrap_loop_state",
        static_cast<std::tuple<int,int,int,int,int> (*)(c_loop_order_periodic&)>(
            __cxxwrap_loop_state
        )
    );

    mod.method(
        "__cxxwrap_restore_loop_state",
        static_cast<void (*)(c_loop_all&,int,int,int,int,int)>(
            __cxxwrap_restore_loop_state
        )
    );
    mod.method(
        "__cxxwrap_restore_loop_state",
        static_cast<void (*)(c_loop_order&,int,int,int,int,int)>(
            __cxxwrap_restore_loop_state
        )
    );
    mod.method(
        "__cxxwrap_restore_loop_state",
        static_cast<void (*)(c_loop_subset&,int,int,int,int,int)>(
            __cxxwrap_restore_loop_state
        )
    );
    mod.method(
        "__cxxwrap_restore_loop_state",
        static_cast<void (*)(c_loop_all_periodic&,int,int,int,int,int)>(
            __cxxwrap_restore_loop_state
        )
    );
    mod.method(
        "__cxxwrap_restore_loop_state",
        static_cast<void (*)(c_loop_order_periodic&,int,int,int,int,int)>(
            __cxxwrap_restore_loop_state
        )
    );
}
