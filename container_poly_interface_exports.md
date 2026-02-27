# Export Status of `container_poly` Class Methods

This table shows which public methods from `container_poly` and its base classes (`container_base`, `wall_list`, `voro_base`, `radius_poly`) are exported to Julia via CxxWrap.jl.

| Method | Signature | Exported As | Status |
|--------|-----------|-------------|--------|
| **Constructor** |
| `container_poly` | `container_poly(double ax_, double bx_, double ay_, double by_, double az_, double bz_, int nx_, int ny_, int nz_, bool xperiodic_, bool yperiodic_, bool zperiodic_, int init_mem)` | `RawContainerPoly(...)` | ✅ Exported |
| `point_inside` | `bool point_inside(double x, double y, double z)` | `__cxxwrap_isinside(con, x, y, z)` | ✅ Exported |
| `draw_domain_gnuplot` | `void draw_domain_gnuplot(FILE *fp = stdout)` | — | ❌ Not exported |
| `draw_domain_gnuplot` | `void draw_domain_gnuplot(const char *filename)` | — | ❌ Not exported |
| `draw_domain_pov` | `void draw_domain_pov(FILE *fp = stdout)` | — | ❌ Not exported |
| `draw_domain_pov` | `void draw_domain_pov(const char *filename)` | — | ❌ Not exported |
| `total_particles` | `int total_particles()` | — | ❌ Not exported |
| `add_wall` | `void add_wall(wall *w)` | — | ❌ Not exported (generic) |
| `add_wall` | `void add_wall(wall &w)` | — | ❌ Not exported (generic) |
| `add_wall` | `void add_wall(wall_list &wl)` | — | ❌ Not exported |
| `point_inside_walls` | `bool point_inside_walls(double x, double y, double z)` | `__cxxwrap_point_inside_walls(con, x, y, z)` | ✅ Exported (lambda) |
| `apply_walls` | `template<class c_class> bool apply_walls(...)` | — | ❌ Not exported |
| **Particle Management** |
| `clear` | `void clear()` | `__cxxwrap_clear!(con)` | ✅ Exported |
| `put` | `void put(int n, double x, double y, double z, double r)` | `__cxxwrap_put!(con, id, x, y, z, r)` | ✅ Exported |
| `put` | `void put(particle_order &vo, int n, double x, double y, double z, double r)` | `__cxxwrap_put!(con, ord, id, x, y, z, r)` | ✅ Exported |
| `import` | `void import(FILE *fp = stdin)` | `import!(con, fptr)` | ✅ Exported |
| `import` | `void import(particle_order &vo, FILE *fp = stdin)` | — | ❌ Not exported |
| `import` | `void import(const char* filename)` | `import!(con, path)` | ✅ Exported |
| `import` | `void import(particle_order &vo, const char *filename)` | — | ❌ Not exported |
| `compute_all_cells` | `void compute_all_cells()` | `compute_all_cells(con)` | ✅ Exported |
| `sum_cell_volumes` | `double sum_cell_volumes()` | `sum_cell_volumes(con)` | ✅ Exported |
| **Cell Computation** |
| `find_voronoi_cell` | `bool find_voronoi_cell(double x, double y, double z, ...)` | `__cxxwrap_find_cell(con, x, y, z)` | ✅ Exported (lambda) |
| `compute_cell` | `template<class v_cell, class c_loop> bool compute_cell(...)` | — | ❌ Not exported |
| `compute_cell` | `template<class v_cell> bool compute_cell(v_cell &c, int ijk, int q)` | — | ❌ Not exported |
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

## Public Data Members

| Member | Type | Exported As | Status |
|--------|------|-------------|--------|
| `max_radius` | `double` | — | ❌ Not exported |

### Notes

- **Template methods** are not exported directly (e.g., `compute_cell`, `apply_walls`, loop-based `draw_*` overloads)
- **`draw_domain_*` methods** are not exported for `container_poly` (only available for `container`)
- **`total_particles`** is not exported for `container_poly` (only available for `container`)
- **`particle_order` overloads** of `import` are not exported
- **`add_wall`** is exported via type-specific lambdas for `wall_sphere`, `wall_cylinder`, `wall_cone`, `wall_plane` only (shared with `container`)
- **`max_radius`** data member is not exported
