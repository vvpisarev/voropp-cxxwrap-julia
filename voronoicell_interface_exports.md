# Export Status of `voronoicell_neighbor` Class Methods

This table shows which public methods from `voronoicell_neighbor` and its base class `voronoicell_base` are exported to Julia via CxxWrap.jl.

| Method | Signature | Exported As | Status |
|--------|-----------|-------------|--------|
| **Constructors/Destructor** |
| `voronoicell_neighbor` | `voronoicell_neighbor()` | `VoronoiCell()` | ✅ Exported |
| `voronoicell_neighbor` | `voronoicell_neighbor(double max_len_sq_)` | `VoronoiCell(max_len_sq)` | ✅ Exported |
| `voronoicell_neighbor` | `template<class c_class> voronoicell_neighbor(c_class &con)` | `VoronoiCell(con)` | ✅ Exported |
| **Assignment** |
| `operator=` | `void operator=(voronoicell_neighbor &c)` | `__cxxwrap_copy!(dest, src)` | ✅ Exported (lambda) |
| **Plane Cutting (with neighbor tracking)** |
| `nplane` | `bool nplane(double x, double y, double z, double rsq, int p_id)` | `__cxxwrap_nplane!(vc, x, y, z, rsq, p_id)` | ✅ Exported |
| `nplane` | `bool nplane(double x, double y, double z, int p_id)` | `__cxxwrap_nplane!(vc, x, y, z, p_id)` | ✅ Exported |
| `plane` | `bool plane(double x, double y, double z, double rsq)` | — | ❌ Not exported |
| `plane` | `bool plane(double x, double y, double z)` | — | ❌ Not exported |
| **Initialization** |
| `init` | `void init(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)` | `__cxxwrap_init!(vc, xmin, xmax, ymin, ymax, zmin, zmax)` | ✅ Exported |
| `init_octahedron` | `void init_octahedron(double l)` | `__cxxwrap_init_octahedron!(vc, l)` | ✅ Exported |
| `init_tetrahedron` | `void init_tetrahedron(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3)` | `__cxxwrap_init_tetrahedron!(vc, ...)` | ✅ Exported |
| **Translation** |
| `translate` | `void translate(double x, double y, double z)` | `__cxxwrap_translate!(vc, x, y, z)` | ✅ Exported |
| **Output (POV-Ray)** |
| `draw_pov` | `void draw_pov(double x, double y, double z, FILE *fp = stdout)` | `__cxxwrap_draw_pov(vc, x, y, z, fptr)` | ✅ Exported |
| `draw_pov` | `void draw_pov(double x, double y, double z, const char *filename)` | `__cxxwrap_draw_pov(vc, x, y, z, filename)` | ✅ Exported |
| `draw_pov_mesh` | `void draw_pov_mesh(double x, double y, double z, FILE *fp = stdout)` | `__cxxwrap_draw_pov_mesh(vc, x, y, z, fptr)` | ✅ Exported |
| `draw_pov_mesh` | `void draw_pov_mesh(double x, double y, double z, const char *filename)` | `__cxxwrap_draw_pov_mesh(vc, x, y, z, filename)` | ✅ Exported |
| **Output (Gnuplot)** |
| `draw_gnuplot` | `void draw_gnuplot(double x, double y, double z, FILE *fp = stdout)` | `__cxxwrap_draw_gnuplot(vc, x, y, z, fptr)` | ✅ Exported |
| `draw_gnuplot` | `void draw_gnuplot(double x, double y, double z, const char *filename)` | `__cxxwrap_draw_gnuplot(vc, x, y, z, filename)` | ✅ Exported |
| **Geometric Properties** |
| `volume` | `double volume()` | `volume(vc)` | ✅ Exported |
| `max_radius_squared` | `double max_radius_squared()` | `max_radius_squared(vc)` | ✅ Exported |
| `total_edge_distance` | `double total_edge_distance()` | `total_edge_distance(vc)` | ✅ Exported |
| `surface_area` | `double surface_area()` | `surface_area(vc)` | ✅ Exported |
| `centroid` | `void centroid(double &cx, double &cy, double &cz)` | `__cxxwrap_centroid(vc)` → `(x, y, z)` | ✅ Exported (lambda) |
| `number_of_faces` | `int number_of_faces()` | `number_of_faces(vc)` | ✅ Exported |
| `number_of_edges` | `int number_of_edges()` | `number_of_edges(vc)` | ✅ Exported |
| **Vertex Information** |
| `vertex_orders` | `void vertex_orders(std::vector<int> &v)` | `__cxxwrap_get_vertex_orders!(v, vc)` | ✅ Exported (lambda) |
| `output_vertex_orders` | `void output_vertex_orders(FILE *fp = stdout)` | — | ❌ Not exported |
| `vertices` | `void vertices(std::vector<double> &v)` | `__cxxwrap_vertices!(v, vc)` | ✅ Exported (lambda) |
| `vertices` | `void vertices(double x, double y, double z, std::vector<double> &v)` | `__cxxwrap_vertices!(v, vc, x, y, z)` | ✅ Exported (lambda) |
| `output_vertices` | `void output_vertices(FILE *fp = stdout)` | — | ❌ Not exported |
| `output_vertices` | `void output_vertices(double x, double y, double z, FILE *fp = stdout)` | — | ❌ Not exported |
| **Face Information** |
| `face_areas` | `void face_areas(std::vector<double> &v)` | `__cxxwrap_face_areas!(v, vc)` | ✅ Exported (lambda) |
| `output_face_areas` | `void output_face_areas(FILE *fp = stdout)` | — | ❌ Not exported |
| `face_orders` | `void face_orders(std::vector<int> &v)` | `__cxxwrap_face_orders!(v, vc)` | ✅ Exported (lambda) |
| `output_face_orders` | `void output_face_orders(FILE *fp = stdout)` | — | ❌ Not exported |
| `face_freq_table` | `void face_freq_table(std::vector<int> &v)` | `__cxxwrap_face_freq_table!(v, vc)` | ✅ Exported (lambda) |
| `output_face_freq_table` | `void output_face_freq_table(FILE *fp = stdout)` | — | ❌ Not exported |
| `face_vertices` | `void face_vertices(std::vector<int> &v)` | `__cxxwrap_face_vertices!(v, vc)` | ✅ Exported (lambda) |
| `output_face_vertices` | `void output_face_vertices(FILE *fp = stdout)` | — | ❌ Not exported |
| `face_perimeters` | `void face_perimeters(std::vector<double> &v)` | `__cxxwrap_face_perimeters!(v, vc)` | ✅ Exported (lambda) |
| `output_face_perimeters` | `void output_face_perimeters(FILE *fp = stdout)` | — | ❌ Not exported |
| **Normal Vectors** |
| `normals` | `void normals(std::vector<double> &v)` | `__cxxwrap_normals!(v, vc)` | ✅ Exported (lambda) |
| `output_normals` | `void output_normals(FILE *fp = stdout)` | — | ❌ Not exported |
| **Custom Output** |
| `output_custom` | `void output_custom(const char *format, FILE *fp = stdout)` | — | ❌ Not exported |
| `output_custom` | `void output_custom(const char *format, int i, double x, double y, double z, double r, FILE *fp)` | — | ❌ Not exported |
| **Plane Cutting** |
| `plane_intersects` | `bool plane_intersects(double x, double y, double z, double rsq)` | `__cxxwrap_plane_intersects(vc, x, y, z, rsq)` | ✅ Exported |
| `plane_intersects_guess` | `bool plane_intersects_guess(double x, double y, double z, double rsq)` | — | ❌ Not exported |
| **Neighbor Information (Virtual)** |
| `neighbors` | `virtual void neighbors(std::vector<int> &v)` | `__cxxwrap_get_neighbors!(v, vc)` | ✅ Exported (lambda) |
| `output_neighbors` | `virtual void output_neighbors(FILE *fp = stdout)` | — | ❌ Not exported |
| **Utilities** |
| `cycle_up` | `inline int cycle_up(int a, int p)` | `__cycle_up(vc, a, p)` | ✅ Exported |
| `cycle_down` | `inline int cycle_down(int a, int p)` | `__cycle_down(vc, a, p)` | ✅ Exported |
| `print_edges` | `void print_edges()` | `print_edges(vc)` | ✅ Exported |

## Public Data Members

| Member | Type | Exported As | Status |
|--------|------|-------------|--------|
| `current_vertices` | `int` | `__get_current_vertices(vc)` | ✅ Exported |
| `p` | `int` | `__get_p(vc)` | ✅ Exported |
| `up` | `int` | `__get_up(vc)` | ✅ Exported |
| `ed` | `int **` | `__cxxwrap_get_ed(vc)` | ✅ Exported |
| `nu` | `int *` | `__cxxwrap_get_nu(vc)` | ✅ Exported |
| `pts` | `double *` | `__cxxwrap_get_pts(vc)` | ✅ Exported |
| `mne` | `int **` | `__cxxwrap_get_mne(vc)` | ✅ Exported |
| `ne` | `int **` | `__cxxwrap_get_ne(vc)` | ✅ Exported |
| `ed[i][j]` | `int` | `__get_ed_ij(vc, i, j)` / `__set_ed_ij!(vc, i, j, k)` | ✅ Exported |
| `tol` | `double` | `__get_tol(vc)` | ✅ Exported |
| `tol_cu` | `double` | — | ❌ Not exported |
| `big_tol` | `double` | — | ❌ Not exported |

### Notes

- **`output_*` methods** taking `FILE*` parameters are generally not exported; use the vector-based alternatives (e.g., `face_areas!`, `vertices!`, `normals!`) instead
- **`plane` methods** (without neighbor tracking) are not exported; use `nplane` with `p_id=0` for similar functionality
- **`plane_intersects_guess`** is not exported; use `plane_intersects` instead
- **Data members** are accessed via getter/setter functions rather than direct access
- **`tol_cu`** and **`big_tol`** data members are not exported
