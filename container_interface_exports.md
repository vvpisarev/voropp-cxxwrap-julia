# Export Status of `container` Class Methods

This table shows which public methods from `container` and its base classes (`container_base`, `wall_list`, `voro_base`) are exported to Julia via CxxWrap.jl.

| Method | Signature | Exported As | Status |
|--------|-----------|-------------|--------|
| **Constructor** |
| `container` | `container(double ax_, double bx_, ...)` | `RawContainer(...)` | ✅ Exported |
| **Inherited from base classes** |
| `point_inside` | `bool point_inside(double x, double y, double z)` | `__cxxwrap_isinside(con, x, y, z)` | ✅ Exported |
| `draw_domain_gnuplot` | `void draw_domain_gnuplot(FILE *fp = stdout)` | — | ❌ Not exported |
| `draw_domain_gnuplot` | `void draw_domain_gnuplot(const char *filename)` | `draw_domain_gnuplot(con, filename)` | ✅ Exported |
| `draw_domain_pov` | `void draw_domain_pov(FILE *fp = stdout)` | — | ❌ Not exported |
| `draw_domain_pov` | `void draw_domain_pov(const char *filename)` | `draw_domain_pov(con, filename)` | ✅ Exported |
| `total_particles` | `int total_particles()` | `total_particles(con)` | ✅ Exported |
| `add_wall` | `void add_wall(wall *w)` | — | ❌ Not exported (generic) |
| `add_wall` | `void add_wall(wall &w)` | — | ❌ Not exported (generic) |
| `add_wall` | `void add_wall(wall_list &wl)` | — | ❌ Not exported |
| `point_inside_walls` | `bool point_inside_walls(double x, double y, double z)` | `__cxxwrap_point_inside_walls(con, x, y, z)` | ✅ Exported |
| `apply_walls` | `template<class c_class> bool apply_walls(...)` | — | ❌ Not exported |
| **Particle Management** |
| `clear` | `void clear()` | `__cxxwrap_clear!(con)` | ✅ Exported |
| `put` | `void put(int n, double x, double y, double z)` | `__cxxwrap_put!(con, id, x, y, z)` | ✅ Exported |
| `put` | `void put(particle_order &vo, int n, double x, double y, double z)` | `__cxxwrap_put!(con, ord, id, x, y, z)` | ✅ Exported |
| `import` | `void import(FILE *fp = stdin)` | `import!(con, fptr)` | ✅ Exported |
| `import` | `void import(particle_order &vo, FILE *fp = stdin)` | — | ❌ Not exported |
| `import` | `void import(const char* filename)` | `import!(con, path)` | ✅ Exported |
| `import` | `void import(particle_order &vo, const char *filename)` | — | ❌ Not exported |
| `compute_all_cells` | `void compute_all_cells()` | `compute_all_cells(con)` | ✅ Exported |
| `sum_cell_volumes` | `double sum_cell_volumes()` | `sum_cell_volumes(con)` | ✅ Exported |
| **Cell Computation** |
| `find_voronoi_cell` | `bool find_voronoi_cell(double x, double y, double z, ...)` | `__cxxwrap_find_cell(con, x, y, z)` | ✅ Exported (lambda) |
| `compute_cell` | `template<class v_cell, class c_loop> bool compute_cell(...)` | `__cxxwrap_compute_cell!(vc, con, cl)` | ✅ Exported (lambda, specific types) |
| **Particle Output (text)** |
| `draw_particles` | `template<class c_loop> void draw_particles(c_loop &vl, FILE *fp)` | — | ❌ Not exported |
| `draw_particles` | `void draw_particles(FILE *fp = stdout)` | `draw_particles(con, fptr)` | ✅ Exported |
| `draw_particles` | `void draw_particles(const char *filename)` | `draw_particles(con, path)` | ✅ Exported |
| **Particle Output (POV-Ray)** |
| `draw_particles_pov` | `template<class c_loop> void draw_particles_pov(c_loop &vl, FILE *fp)` | — | ❌ Not exported |
| `draw_particles_pov` | `void draw_particles_pov(FILE *fp = stdout)` | `draw_particles_pov(con, fptr)` | ✅ Exported |
| `draw_particles_pov` | `void draw_particles_pov(const char *filename)` | `draw_particles_pov(con, path)` | ✅ Exported |
| **Cell Output (Gnuplot)** |
| `draw_cells_gnuplot` | `template<class c_loop> void draw_cells_gnuplot(c_loop &vl, FILE *fp)` | — | ❌ Not exported |
| `draw_cells_gnuplot` | `void draw_cells_gnuplot(FILE *fp = stdout)` | `draw_cells_gnuplot(con, fptr)` | ✅ Exported |
| `draw_cells_gnuplot` | `void draw_cells_gnuplot(const char *filename)` | `draw_cells_gnuplot(con, path)` | ✅ Exported |
| **Cell Output (POV-Ray)** |
| `draw_cells_pov` | `template<class c_loop> void draw_cells_pov(c_loop &vl, FILE *fp)` | — | ❌ Not exported |
| `draw_cells_pov` | `void draw_cells_pov(FILE *fp = stdout)` | `draw_cells_pov(con, fptr)` | ✅ Exported |
| `draw_cells_pov` | `void draw_cells_pov(const char *filename)` | `draw_cells_pov(con, path)` | ✅ Exported |

### Notes

- **Template methods** are generally not exported directly; some have specialized lambda wrappers for specific types (e.g., `compute_cell`, `find_voronoi_cell`)
- **`add_wall`** is exported via type-specific lambdas for `wall_sphere`, `wall_cylinder`, `wall_cone`, `wall_plane` only
- **`FILE*` overloads** for `draw_domain_*` methods are not exported (only `const char*` versions)
- **`particle_order` overloads** of `import` are not exported
- **Loop-based overloads** (`c_loop&` parameter) are not exported
