#include "scene.hpp"
#define JC_VORONOI_IMPLEMENTATION
#include "jc_voronoi.h"
#define JC_VORONOI_CLIP_IMPLEMENTATION
#include "jc_voronoi_clip.h"

const float SNOW_HEIGHT = 0.8f;
const float SNOW_STD = 0.2f;

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
    std::default_random_engine e(time(NULL));

	create_voronoi(2000);

	create_ocean_border();

	compute_heights();

	compute_waterdists();

	// Compute the elevation of the centers of the polygons
	smooth_centers();

	/// Laplacian smoothing
	/*for (int i = 0; i < 5; i++) {
		laplacian_smoothing();
	}*/

	add_snow_biotope();

	// Create a visual frame representing the coordinate system
	global_frame.initialize(mesh_primitive_frame(), "Frame");
	environment.camera.axis = camera_spherical_coordinates_axis::z;
	environment.camera.look_at({ 5.0f,5.0f,15.0f }, { 5,5,0 });

	// Create a mesh for the centers
	particle_sphere.initialize(mesh_primitive_sphere(0.05f));
	particle_sphere.shading.color = { 1,0,0 };

	create_terrain();

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
	draw(sea, environment);
	draw(terrain, environment);
}

/// Displays the GUI elements.
void scene_structure::display_gui() {
	ImGui::Checkbox("Frame", &gui.display_frame);
}

/// Creates a terrain mesh based on the Voronoi diagram.
void scene_structure::create_terrain() {
	mesh sea_mesh = mesh_primitive_grid({-20,-20,-0.01f},{20,-20,-0.01f},{20,20,-0.01f},{-20,20,-0.01f});
	// Initialize and sets the right color for the terrain
	sea.initialize(sea_mesh, "sea");
	sea.shading.color = { 0,0.412f,0.58f };
	sea.shading.phong.specular = 0;

    // Terrain geometry
	mesh terrain_mesh;
    for(int k = 0; k < N; k++) {
		// The following convention is being used:
		// the first point is the centers of the polygon,
		// then the following points are its corners.
		mesh polygon_mesh;
		const int center_idx = 0;

		vec3 color;
		if (biotopes[k] == Biotope::Ocean)
			color = vec3(0,0.412f,0.58f);
		else if (biotopes[k] == Biotope::Lake)
			color = vec3(0.2f,0.59f,0.885f);
		else if (biotopes[k] == Biotope::Snow)
			color = vec3(1.0f,1.0f,1.0f);
		else
			color = vec3(0.6f,0.85f,0.5f);

		polygon_mesh.position.push_back(centers[k]);
		polygon_mesh.color.push_back(color);
		for(int& corner: mycorners[k]) {
			polygon_mesh.position.push_back(corners[corner]);
			polygon_mesh.color.push_back(color);
		}

		// Generate triangle organization
		// The triangles are constructed as follows: a triangle is made up of
		// the center of a polygon, and the vectices of an edge of this polygon
		for(int i = 0; i < mycorners[k].size(); i++) {
			int corner_1 = 1+(i%mycorners[k].size());
			int corner_2 = 1+((i+1)%mycorners[k].size());
        	uint3 triangle = {center_idx, corner_1, corner_2};
        	polygon_mesh.connectivity.push_back(triangle);
		}

		// Apply defaults
		polygon_mesh.fill_empty_field();
		terrain_mesh.push_back(polygon_mesh);
	}

	// Initialize and sets the right color for the land
	terrain.initialize(terrain_mesh, "Terrain");
	terrain.shading.phong.specular = 0.0f; // non-specular land material
}

/// Computes the height, and find the connected component
/// of the oceans at the same time.
void scene_structure::compute_heights() {
	/// Find the source of the shortest path algorithm.
	int source = 0;
	for (int i = 0; i < N_corners; i++) {
		corners[i].z = INFINITY;
		// If there is a ocean polygon amongst all the polygon
		// this corners touches, sets the height to 0 
		for (int& poly_idx: touches[i]) {
			if (biotopes[poly_idx] == Biotope::Ocean) {
				corners[i].z = 0;
			}
		}
	}
	// In case of, make sure that the source is at the sea level
	corners[source].z = 0;

	// Build the priority queue
	std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::function<bool(std::pair<int,int>, std::pair<int,int>)>> Q(CompareHeights);
	for (int v = 0; v < N_corners; v++) {
		Q.push(std::pair<int,int>(v, corners[v].z));
	}
	
	while (!Q.empty()) {
		int u = Q.top().first;
		Q.pop();
		for (auto& adjacent: adjacents[u]) {
			int v = std::get<0>(adjacent);
			int poly_1 = std::get<1>(adjacent);
			int poly_2 = std::get<1>(adjacent);
			float dist;
			if (biotopes[poly_1] == Biotope::Ocean) {
				if (biotopes[poly_2] == Biotope::Lake) biotopes[poly_2] = Biotope::Ocean;
				dist = 0;
			} else if (biotopes[poly_2] == Biotope::Ocean) {
				if (biotopes[poly_1] == Biotope::Lake) biotopes[poly_1] = Biotope::Ocean;
				dist = 0;
			} else if (biotopes[poly_1] == Biotope::Lake || biotopes[poly_2] == Biotope::Lake) {
				dist = 0;
			} else {
				vec2 corner_u = { corners[u].x, corners[u].y };
				vec2 corner_v = { corners[v].x, corners[v].y };
				dist = norm(corner_u - corner_v);
			}

			float alt = corners[u].z + dist;
			if (alt < corners[v].z) {
				corners[v].z = alt;
				Q.push(std::pair<int,int>(v, corners[v].z));
			}
		}
	}
	std::cout << "Heights of corners computed";
}

/// Computes the smallest distance to a water source.
void scene_structure::compute_waterdists() {
	waterdists.resize(N_corners);
	// Compute the source of the algorithm
	// and initializes the distances
	int source = 0;
	for (int i = 0; i < N_corners; i++) {
		waterdists[i] = INFINITY;
		// If there is a ocean polygon amongst all the polygon
		// this corners touches, sets the height to 0 
		for (int& poly_idx: touches[i]) {
			if (biotopes[poly_idx] == Biotope::Ocean || biotopes[poly_idx] == Biotope::Lake) {
				waterdists[i] = 0;
			}
		}
	}
	// In case of, make sure that the source is at the sea level
	waterdists[source] = 0;

	// Build the priority queue
	std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::function<bool(std::pair<int,int>, std::pair<int,int>)>> Q(CompareHeights);
	for (int v = 0; v < N_corners; v++) {
		Q.push(std::pair<int,int>(v, waterdists[v]));
	}
	
	while (!Q.empty()) {
		int u = Q.top().first;
		Q.pop();
		for (auto& adjacent: adjacents[u]) {
			int v = std::get<0>(adjacent);
			int poly_1 = std::get<1>(adjacent);
			int poly_2 = std::get<1>(adjacent);
			float alt;
			if (biotopes[poly_1] == Biotope::Ocean || biotopes[poly_1] == Biotope::Lake ||
				biotopes[poly_2] == Biotope::Ocean || biotopes[poly_2] == Biotope::Lake) {
				alt = 0;
			} else {
				vec2 corner_u = { corners[u].x, corners[u].y };
				vec2 corner_v = { corners[v].x, corners[v].y };
				alt = waterdists[u] + norm(corner_u - corner_v);
			}

			if (alt < waterdists[v]) {
				waterdists[v] = alt;
				Q.push(std::pair<int,int>(v, waterdists[v]));
			}
		}
	}

	std::cout << "Distance to water computed";
}

/// Creates a new Voronoi structure with `n` clusters
/// and updates the class fields accordingly.
void scene_structure::create_voronoi(int n) {
	N = n;
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

	// Always free up memory at the end
    jcv_diagram_free( &diagram );
}

/// Creates an ocean border around the island.
void scene_structure::create_ocean_border() {
	for(int idx = 0; idx < N; idx++) {
		for (auto& corner_idx: mycorners[idx]) {
			vec3& corner = corners[corner_idx];
			if (corner.x <= 0.2f || corner.y <= 0.2f || corner.x >= 9.8f || corner.y  >= 9.8f) {
				  biotopes[idx] = Biotope::Ocean;
			}
		}
		// Old values:
		// * threshold: 0.5f
		// * noise_perlin: 6, 0.37f, 2.2f
		// * noise < threshold
		float threshold = 1.0f;
		if (noise_perlin({centers[idx].x, centers[idx].y}, 2, 0.1f, 1.0f) > threshold) {
			biotopes[idx] = Biotope::Lake;
		}
	}
}

/// Places the center of the polygons at the mean of the elevation its corners.
void scene_structure::smooth_centers() {
	for (int v = 0; v < N; v++) {
		float mean = 0.0f;
		int n = 0;
		for (int& corner_idx: mycorners[v]) {
			mean += corners[corner_idx].z;
			n++;
		}
		centers[v].z = mean/((float) n);
	}
}

/// Applies Laplacian smoothing to both the centers and the corners of each polygon.
void scene_structure::laplacian_smoothing() {
	for (int v = 0; v < N_corners; v++) {
		float mean = 0.0f;
		int n = 0;
		for (auto& adjacent: adjacents[v]) {
			mean += corners[std::get<0>(adjacent)].z;
			n++;
		}
		for (int& poly: touches[v]) {
			mean += centers[poly].z;
			n++;
		}
		corners[v].z = mean/((float) n);
	}

	smooth_centers();
}

/// Adds all the other biotopes.
void scene_structure::add_biotopes() {
	// We define snow biotopes based on elevation
    std::normal_distribution<double> distSnow(SNOW_HEIGHT,SNOW_STD); 
	for (int idx = 0; idx < N; idx++) {
		if (centers[idx].z > distSnow(e) && biotopes[idx] == Biotope::Land) {
			biotopes[idx] = Biotope::Snow;
		}
	}	
}