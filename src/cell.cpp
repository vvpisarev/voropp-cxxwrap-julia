void export_voronoicell_methods(jlcxx::Module& mod)
{
    using namespace voro;

    // voronoicell_base methods
    mod.method(
        "__cxxwrap_current_vertices",
        static_cast<int (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.current_vertices; }
        )
    );
    mod.method(
        "__cxxwrap_p",
        static_cast<int (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.p; }
        )
    );
    mod.method(
        "__cxxwrap_up",
        static_cast<int (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.up; }
        )
    );
    mod.method(
        "__cxxwrap_ed",
        static_cast<int** (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.ed; }
        )
    );
    mod.method(
        "__cxxwrap_nu",
        static_cast<int* (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.nu; }
        )
    );
    mod.method(
        "__cxxwrap_pts",
        static_cast<double* (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.pts; }
        )
    );
    mod.method(
        "__cxxwrap_tol",
        static_cast<double (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.tol; }
        )
    );
    mod.method(
        "__cxxwrap_translate!",
        static_cast<void (*)(voronoicell_neighbor&, double, double, double)>(
            [] (voronoicell_neighbor& v, double x, double y, double z)
            {
                v.translate(x, y, z);
            }
        )
    );
    mod.method(
        "__cxxwrap_draw_pov",
        static_cast<void (*)(FILE*, voronoicell_neighbor&, double, double, double)>(
            [] (FILE* fp, voronoicell_neighbor& vc, double x, double y, double z) {
                vc.draw_pov(x, y, z, fp);
            }
        )
    );
    mod.method(
        "__cxxwrap_draw_pov",
        static_cast<void (*)(const char*, voronoicell_neighbor&, double, double, double)>(
            [] (const char* fp, voronoicell_neighbor& vc, double x, double y, double z) {
                vc.draw_pov(x, y, z, fp);
            }
        )
    );
    mod.method(
        "__cxxwrap_draw_pov_mesh",
        static_cast<void (*)(FILE*, voronoicell_neighbor&, double, double, double)>(
            [] (FILE* fp, voronoicell_neighbor& vc, double x, double y, double z) {
                vc.draw_pov_mesh(x, y, z, fp);
            }
        )
    );
    mod.method(
        "__cxxwrap_draw_pov_mesh",
        static_cast<void (*)(const char*, voronoicell_neighbor&, double, double, double)>(
            [] (const char* fp, voronoicell_neighbor& vc, double x, double y, double z) {
                vc.draw_pov_mesh(x, y, z, fp);
            }
        )
    );
    mod.method(
        "__cxxwrap_draw_gnuplot",
        static_cast<void (*)(FILE*, voronoicell_neighbor&, double, double, double)>(
            [] (FILE* fp, voronoicell_neighbor& vc, double x, double y, double z) {
                vc.draw_gnuplot(x, y, z, fp);
            }
        )
    );
    mod.method(
        "__cxxwrap_draw_gnuplot",
        static_cast<void (*)(const char*, voronoicell_neighbor&, double, double, double)>(
            [] (const char* fp, voronoicell_neighbor& vc, double x, double y, double z) {
                vc.draw_gnuplot(x, y, z, fp);
            }
        )
    );
    mod.method(
        "__cxxwrap_volume",
        static_cast<double (*)(voronoicell_neighbor&)>(
            [](voronoicell_neighbor& vc) { return vc.volume(); }
        )
    );
    mod.method(
        "__cxxwrap_max_radius_squared",
        static_cast<double (*)(voronoicell_neighbor&)>(
            [](voronoicell_neighbor& vc) { return vc.max_radius_squared(); }
        )
    );
    mod.method(
        "__cxxwrap_total_edge_distance",
        static_cast<double (*)(voronoicell_neighbor&)>(
            [](voronoicell_neighbor& vc) { return vc.total_edge_distance(); }
        )
    );
    mod.method(
        "__cxxwrap_surface_area",
        static_cast<double (*)(voronoicell_neighbor&)>(
            [](voronoicell_neighbor& vc) { return vc.surface_area(); }
        )
    );
    mod.method(
        "__cxxwrap_centroid",
        static_cast<vec3d (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v)
            {
                vec3d rc;
                v.centroid(rc.x, rc.y, rc.z);
                return rc;
            }
        )
    );
    mod.method(
        "__cxxwrap_number_of_faces",
        static_cast<int (*)(voronoicell_neighbor&)>(
            [](voronoicell_neighbor& v) { return v.number_of_faces(); }
        )
    );
    mod.method(
        "__cxxwrap_number_of_edges",
        static_cast<int (*)(voronoicell_neighbor&)>(
            [](voronoicell_neighbor& v) { return v.number_of_edges(); }
        )
    );
    mod.method(
        "__cxxwrap_vertex_orders!",
        static_cast<void (*)(std::vector<int>&, voronoicell_neighbor&)>(
            [] (std::vector<int>& v, voronoicell_neighbor& vc)
            {
                vc.vertex_orders(v);
            }
        )
    );
    mod.method(
        "__cxxwrap_vertices!",
        static_cast<void (*)(std::vector<double>&, voronoicell_neighbor&)>(
            [] (std::vector<double>& v, voronoicell_neighbor& vc)
            {
                vc.vertices(v);
            }
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
        "__cxxwrap_face_areas!",
        static_cast<void (*)(std::vector<double>&, voronoicell_neighbor&)>(
            [] (std::vector<double>& v, voronoicell_neighbor& vc)
            {
                vc.face_areas(v);
            }
        )
    );
    mod.method(
        "__cxxwrap_face_orders!",
        static_cast<void (*)(std::vector<int>&, voronoicell_neighbor&)>(
            [] (std::vector<int>& v, voronoicell_neighbor& vc)
            {
                vc.face_orders(v);
            }
        )
    );
    mod.method(
        "__cxxwrap_face_freq_table!",
        static_cast<void (*)(std::vector<int>&, voronoicell_neighbor&)>(
            [] (std::vector<int>& v, voronoicell_neighbor& vc)
            {
                vc.face_freq_table(v);
            }
        )
    );
    mod.method(
        "__cxxwrap_face_vertices!",
        static_cast<void (*)(std::vector<int>&, voronoicell_neighbor&)>(
            [] (std::vector<int>& v, voronoicell_neighbor& vc)
            {
                vc.face_vertices(v);
            }
        )
    );
    mod.method(
        "__cxxwrap_face_perimeters!",
        static_cast<void (*)(std::vector<double>&, voronoicell_neighbor&)>(
            [] (std::vector<double>& v, voronoicell_neighbor& vc)
            {
                vc.face_perimeters(v);
            }
        )
    );
    mod.method(
        "__cxxwrap_normals!",
        static_cast<void (*)(std::vector<double>&, voronoicell_neighbor&)>(
            [] (std::vector<double>& v, voronoicell_neighbor& vc)
            {
                vc.normals(v);
            }
        )
    );
    mod.method(
        "__cxxwrap_cycle_up",
        static_cast<int (*)(voronoicell_neighbor&, int, int)>(
            [](voronoicell_neighbor& vc, int a, int p) { return vc.cycle_up(a, p); }
        )
    );
    mod.method(
        "__cxxwrap_cycle_down",
        static_cast<int (*)(voronoicell_neighbor&, int, int)>(
            [](voronoicell_neighbor& vc, int a, int p) { return vc.cycle_down(a, p); }
        )
    );
    // voronoicell_neighbor methods
    mod.method(
        "__cxxwrap_mne",
        static_cast<int** (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.mne; }
        )
    );
    mod.method(
        "__cxxwrap_ne",
        static_cast<int** (*)(voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& v) { return v.ne; }
        )
    );
    mod.method(
        "__cxxwrap_copyto!",
        static_cast<void (*)(voronoicell_neighbor&, voronoicell_neighbor&)>(
            [] (voronoicell_neighbor& dest, voronoicell_neighbor& src) {
                dest = src;
            }
        )
    );
    mod.method(
        "__cxxwrap_nplane!",
        static_cast<bool (*)(voronoicell_neighbor&, double, double, double, double, int)>(
            [] (voronoicell_neighbor& vc, double x,double y,double z,double rsq,int p_id) 
            {
                return vc.nplane(x, y, z, rsq, p_id);
            }
        )
    );
    mod.method(
        "__cxxwrap_nplane!",
        static_cast<bool (*)(voronoicell_neighbor&, double, double, double, int)>(
            [] (voronoicell_neighbor &vc, double x,double y,double z,int p_id) 
            {
                return vc.nplane(x, y, z, p_id);
            }
        )
    );
    mod.method(
        "__cxxwrap_plane!",
        static_cast<bool (*)(voronoicell_neighbor&, double, double, double, double)>(
            [] (voronoicell_neighbor &vc, double x,double y,double z,double rsq) 
            {
                return vc.plane(x, y, z, rsq);
            }
        )
    );
    mod.method(
        "__cxxwrap_plane!",
        static_cast<bool (*)(voronoicell_neighbor&, double, double, double)>(
            [] (voronoicell_neighbor &vc, double x,double y,double z) 
            {
                return vc.plane(x, y, z);
            }
        )
    );
    mod.method(
        "__cxxwrap_init!",
        static_cast<void (*)(voronoicell_neighbor&, double, double, double, double, double, double)>(
            [] (voronoicell_neighbor& vc, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
            {
                vc.init(xmin, xmax, ymin, ymax, zmin, zmax);
            }
        )
    );
    mod.method(
        "__cxxwrap_init_octahedron!",
        static_cast<void (*)(voronoicell_neighbor&, double)>(
            [] (voronoicell_neighbor& vc, double l)
            {
                vc.init_octahedron(l);
            }
        )
    );
    mod.method(
        "__cxxwrap_init_tetrahedron!",
        static_cast<void (*)(voronoicell_neighbor&, double, double, double, double, double, double, double, double, double, double, double, double)>(
            [] (voronoicell_neighbor& vc, double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)
            {
                vc.init_tetrahedron(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);
            }
        )
    );
    mod.method(
        "__cxxwrap_neighbors!",
        static_cast<void (*)(std::vector<int>&, voronoicell_neighbor&)>(
            [] (std::vector<int> &v, voronoicell_neighbor &vc)
            {
                vc.neighbors(v);
            }
        )
    );
}
