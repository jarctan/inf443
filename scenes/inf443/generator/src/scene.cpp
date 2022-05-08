#include "scene.hpp"
#define JC_VORONOI_IMPLEMENTATION
#include "jc_voronoi.h"
#define JC_VORONOI_CLIP_IMPLEMENTATION
#include "jc_voronoi_clip.h"

using namespace cgp;

/// Place the points at the centers (centroid) of their respective clusters
/// in the Voronoi diagram.
void relax_points(const jcv_diagram* diagram, jcv_point* points) {
    const jcv_site* sites = jcv_diagram_get_sites(diagram);
    for( int i = 0; i < diagram->numsites; ++i )
    {
        const jcv_site* site = &sites[i];
        jcv_point sum = {0, 0};
        int count = 0;

        const jcv_graphedge* edge = site->edges;

        while( edge )
        {
            sum.x += edge->pos[0].x;
            sum.y += edge->pos[0].y;
            ++count;
            edge = edge->next;
        }

        points[site->index].x = sum.x / count;
        points[site->index].y = sum.y / count;
    }
}

/// This function compares the heights of two pairs of (_,height)
/// and returns the true if the first one is higher than the second one.
/// This is only used to create the priority queue.
bool scene_structure::CompareHeights(const std::pair<int,int> &a, const std::pair<int,int> &b) {
    return a.second > b.second;
}

/// Find a point in an array, and return the index of this point.
int find_pt(std::vector<vec3> l, vec3 el) {
	for(int i = 0; i < l.size(); i++) {
		if (norm(el - l[i]) <= 0.0001f) {
			return i;
		}
	}
	return -1;
}

/// This function is called only once at the beginning of the program
/// and initializes the meshes, diagrams and other structures.
void scene_structure::initialize() {
	// Number of clusters
	N = 3000;
	centers.resize(N);
	biotopes.resize(N, Biotope::Land);
	neighbors.resize(N);
	mycorners.resize(N);

	jcv_diagram diagram;
    memset(&diagram, 0, sizeof(jcv_diagram));
	jcv_point* points = (jcv_point*) malloc(N*sizeof(jcv_point));
	for (int i = 0; i < N; i++) {
		jcv_point pt;
		pt.x = ((float) rand()) * 10.0f / RAND_MAX;
		pt.y = (float) rand() * 10.0f / RAND_MAX;
		points[i] = pt;
	}

    // Create a clipper, which is a rectangle
	jcv_clipping_polygon polygon;
    polygon.num_points = 4;
    polygon.points = (jcv_point*)malloc(sizeof(jcv_point)*(size_t)polygon.num_points);

    polygon.points[0].x = 0;
    polygon.points[0].y = 0;
    polygon.points[1].x = 10;
    polygon.points[1].y = 0;
    polygon.points[2].x = 10;
    polygon.points[2].y = 10;
    polygon.points[3].x = 0;
    polygon.points[3].y = 10;

    jcv_clipper polygonclipper;
    polygonclipper.test_fn = jcv_clip_polygon_test_point;
    polygonclipper.clip_fn = jcv_clip_polygon_clip_edge;
    polygonclipper.fill_fn = jcv_clip_polygon_fill_gaps;
    polygonclipper.ctx = &polygon;

	// Run n iterations of Voronoi (find Voronoi clusters, then place points at centroids, repeat)
	for (int i = 0; i < 1; i++) {
		jcv_diagram_generate(N, points, 0, &polygonclipper, &diagram );
		relax_points(&diagram, points);
	}

	segments_drawable::default_shader = curve_drawable::default_shader;
	segment.initialize({ {0,0,0},{1,0,0} });

	// Store the centers
	for(int i = 0; i < N; i++) {
		centers[i] = {points[i].x, points[i].y, 0};
	}

	// Store the corners and neighbors
    const jcv_site* sites = jcv_diagram_get_sites(&diagram);
    for(int i = 0; i < diagram.numsites; i++) {
        const jcv_site* site = &sites[i];
        const jcv_graphedge* edge = site->edges;
		int idx = site->index;
		std::vector<std::tuple<int,int,int>> l;

        while( edge ) {
			vec3 v = {edge->pos[0].x, edge->pos[0].y, 0.0};

			// For each point in the edge, retrieve its index in the corners list.
			// If it does not already exist, add it to the list (the index will then
			// be the size of the list).
			int pt1 = find_pt(corners, v);
			if (pt1 == -1) {
				pt1 = corners.size();
				corners.push_back(v);
				touches.push_back({});
				adjacents.push_back({});
			}
			// If there's some logic here, the list of corners is exactly
			// the list of the first vertex of each edge in the polygon
			mycorners[idx].push_back(pt1);

			v = {edge->pos[1].x, edge->pos[1].y, 0.0};
			int pt2 = find_pt(corners, v);
			if (pt2 == -1) {
				pt2 = corners.size();
				corners.push_back(v);
				touches.push_back({});
				adjacents.push_back({});
			}

			int idx_neigh = edge->neighbor ? edge->neighbor->index : -1;
			adjacents[pt1].push_back(std::make_tuple(pt2, idx, idx_neigh));
			adjacents[pt2].push_back(std::make_tuple(pt1, idx, idx_neigh));

			touches[pt1].push_back(idx);
			touches[pt2].push_back(idx);

			// Add the edge in the form of (neighbor, idx of 1st corner, idx of 2nd corner)
			if (edge->neighbor)
				l.push_back(std::make_tuple(edge->neighbor->index, pt1, pt2));
				
            edge = edge->next;
		}

		neighbors[idx] = l;
	}
	N_corners = corners.size();

	std::cout << "End of preprocessing Voronoi diagrams";

	// Create the biotopes
	for(int i = 0; i < diagram.numsites; i++) {
        const jcv_site* site = &sites[i];
        const jcv_graphedge* edge = site->edges;
		int idx = site->index;
		
		while (edge) {
			vec3 p1 = {edge->pos[0].x, edge->pos[0].y, 0 };
			vec3 p2 = {edge->pos[1].x, edge->pos[1].y, 0 };
			if (p1.x <= 0.2f || p1.y <= 0.2f || p2.x <= 0.2f || p2.y <= 0.2f
			  ||p1.x >= 9.8f || p2.x >= 9.8f || p2.x >= 9.8f || p2.y >= 9.8f) {
				  biotopes[idx] = Biotope::Ocean;
			}
				
            edge = edge->next;
		}
		//Old values:
		// * threshold: 0.5f
		// * noise_perlin: 6, 0.37f, 2.2f
		// * noise < threshold
		float threshold = 1.0f;
		if (noise_perlin({centers[idx].x, centers[idx].y}, 6, 0.7f, 0.6f) < threshold) {
			biotopes[idx] = Biotope::Ocean;
		}
	}
	
	// Compute the height
	heights.resize(N_corners);
	// Find the source of the shortest path algorithm
	int source = 0;
	for (int i = 0; i < N_corners; i++) {
		heights[i] = INFINITY;
		// If there is a ocean polygon amongst all the polygon
		// this corners touches, sets the height to 0 
		for (int& poly_idx: touches[i]) {
			if (biotopes[poly_idx] == Biotope::Ocean) {
				heights[i] = 0;
			}
		}
	}
	// In case of, make sure that the source is at the sea level
	heights[source] = 0;

	// Build the priority queue
	std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::function<bool(std::pair<int,int>, std::pair<int,int>)>> Q(CompareHeights);
	for (int v = 0; v < N_corners; v++) {
		Q.push(std::pair<int,int>(v, heights[v]));
	}
	
	while (!Q.empty()) {
		int u = Q.top().first;
		Q.pop();
		for (auto& adjacent: adjacents[u]) {
			int v = std::get<0>(adjacent);
			float dist = norm(corners[u] - corners[v]);
			float alt = heights[u] + dist;
			if (alt < heights[v]) {
				heights[v] = alt;
				Q.push(std::pair<int,int>(v, heights[v]));
			}
		}
	}

	/// Corners heights have just been computed
	for (int v = 0; v < N_corners; v++) {
		corners[v].z = heights[v];
	}

	// We now compute the elevation of the centers
	// of the polygons
	for (int v = 0; v < N; v++) {
		// The elevation of the center is
		// the mean of the elevation of its corners
		float mean = 0.0f;
		int n = 0;
		for (int& corner_idx: mycorners[v]) {
			mean += corners[corner_idx].z;
			n++;
		}
		centers[v].z = mean/((float) n);
	}

	// Create a visual frame representing the coordinate system
	global_frame.initialize(mesh_primitive_frame(), "Frame");
	environment.camera.axis = camera_spherical_coordinates_axis::z;
	environment.camera.look_at({ 5.0f,5.0f,10.0f }, { 5,5,0 });

	// Create a mesh for the centers
	particle_sphere.initialize(mesh_primitive_sphere(0.05f));
	particle_sphere.shading.color = { 1,0,0 };

	// Create the terrain
	create_terrain();

	// Always free up memory at the end
    jcv_diagram_free( &diagram );

	// Timer scale
	timer.scale = 0.5f;
}

/// Draws (and displays) a segment between a and b.
void scene_structure::draw_segment(vec3 const& a, vec3 const& b) {
	segment.update({ a, b });
	draw(segment, environment);
}

// This function is constantly called at every frame
void scene_structure::display() {
	// Set the light to the current position of the camera
	environment.light = environment.camera.position();

	// conditional display of the global frame (set via the GUI)
	if (gui.display_frame)
		draw(global_frame, environment);

	// Update the current time
	timer.update();

	// Draw the lines between the clusters, and the borders
	draw_segment({0,0,0},{0,10,0});
	draw_segment({0,10,0},{10,10,0});
	draw_segment({10,10,0},{10,0,0});
	draw_segment({10,0,0},{0,0,0});
	draw(terrain, environment);

	if (gui.display_wireframe);
}

void scene_structure::display_gui() {
	ImGui::Checkbox("Frame", &gui.display_frame);
	ImGui::Checkbox("Wireframe", &gui.display_wireframe);
}

void scene_structure::create_terrain()
{
	// The following convention is being used:
	// the first N points are the centers of the polygons,
	// then the N_corners following are the points of the corners.
    mesh terrain_mesh;
    terrain_mesh.position.resize(N + N_corners);

    // Fill terrain geometry
    for(int k = 0; k < N; k++)
		terrain_mesh.position[k] = centers[k];
    for(int k = 0; k < N_corners; k++)
		terrain_mesh.position[N+k] = corners[k];

    // Generate triangle organization
	// The triangles are constructed as follows: a triangle is made up of
	// the center of a polygon, and the vectices of an edge of this polygon
    for(int k = 0; k < N; k++) {
		for (auto& neighbor: neighbors[k]) {
        	uint3 triangle = {k, N+std::get<1>(neighbor), N+std::get<2>(neighbor)};
        	terrain_mesh.connectivity.push_back(triangle);
		}
    }

	// Apply defaults
	terrain_mesh.fill_empty_field(); 

	// Initialize and sets the right color for the terrain
	terrain.initialize(terrain_mesh, "terrain");
	terrain.shading.color = { 0.6f,0.85f,0.5f };
	terrain.shading.phong.specular = 0.0f; // non-specular terrain material
}