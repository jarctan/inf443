#include "scene.hpp"

const float SNOW_HEIGHT = 0.8f;
const float SNOW_STD = 0.2f;

using namespace cgp;
using namespace std::chrono;

/// This function is called only once at the beginning of the program
/// and initializes the meshes, diagrams and other structures.
void scene_structure::initialize() {
	randeng = std::default_random_engine(time(NULL));

	auto start = high_resolution_clock::now();
	create_voronoi(6000);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << "Voronoi diagram in " << duration.count() << "ms [OK]" << std::endl;

	start = high_resolution_clock::now();
	create_ocean_border();
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	std::cout << "Ocean border in " << duration.count() << "ms [OK]" << std::endl;

	start = high_resolution_clock::now();
	compute_heights();
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	std::cout << "Compute height in " << duration.count() << "ms [OK]" << std::endl;

	start = high_resolution_clock::now();
	compute_waterdists();
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	std::cout << "Water distance in " << duration.count() << "ms [OK]" << std::endl;

	// Compute the elevation of the centers of the polygons
	start = high_resolution_clock::now();
	smooth_centers();
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	std::cout << "Smooth centers in " << duration.count() << "ms [OK]" << std::endl;

	start = high_resolution_clock::now();
	add_biotopes();
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	std::cout << "Biotopes in " << duration.count() << "ms [OK]" << std::endl;

	// Create a visual frame representing the coordinate system
	global_frame.initialize(mesh_primitive_frame(), "Frame");
	environment.camera.axis = camera_spherical_coordinates_axis::z;
	environment.camera.look_at({ 5.0f,5.0f,15.0f }, { 5,5,0 });;

	start = high_resolution_clock::now();
	create_terrain();
	stop = high_resolution_clock::now();
	std::cout << "Terrain in " << duration.count() << "ms [OK]" << std::endl;

	timer.scale = 0.5f;
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
    for (int k = 0; k < N; k++) {
		// The following convention is being used:
		// the first point is the centers of the polygon,
		// then the following points are its corners.
		mesh polygon_mesh;
		const int center_idx = 0;

		vec4 color;
		if (biotopes[k] == Biotope::Ocean)
			color = vec4(0,0.412f,0.58f,0.0f);
		else if (biotopes[k] == Biotope::Lake)
			color = vec4(0.2f,0.59f,0.885f,1.0f);
		else if (biotopes[k] == Biotope::Snow)
			color = vec4(1.0f,1.0f,1.0f,1.0f);
		else if (biotopes[k] == Biotope::Tundra)
			color = vec4(0.867f,0.867f,0.733f,1.0f);
		else if (biotopes[k] == Biotope::Bare)
			color = vec4(0.733f,0.733f,0.733f,1.0f);
		else if (biotopes[k] == Biotope::Scorched)
			color = vec4(0.6f,0.6f,0.6f,1.0f);
		else if (biotopes[k] == Biotope::Taiga)
			color = vec4(0.8f,0.831f,0.733f,1.0f);
		else if (biotopes[k] == Biotope::Shrubland)
			color = vec4(0.769f,0.8f,0.733f,1.0f);
		else if (biotopes[k] == Biotope::TempDesert)
			color = vec4(0.894f,0.91f,0.792f,1.0f);
		else if (biotopes[k] == Biotope::TempRainForest)
			color = vec4(0.643f,0.769f,0.659f,1.0f);
		else if (biotopes[k] == Biotope::TempDeciduousForest)
			color = vec4(0.706f,0.788f,0.663f,1.0f);
		else if (biotopes[k] == Biotope::Grassland)
			color = vec4(0.769f,0.831f,0.667f,1.0f);
		else if (biotopes[k] == Biotope::TropicalRainForest)
			color = vec4(0.612f,0.733f,0.663f,1.0f);
		else if (biotopes[k] == Biotope::TropicalSeasonalForest)
			color = vec4(0.663f,0.8f,0.643f,1.0f);
		else if (biotopes[k] == Biotope::SubtropicalDesert)
			color = vec4(0.914f,0.867f,0.78f,1.0f);
		else
			color = vec4(0.6f,0.85f,0.5f,0.0f);

		polygon_mesh.position.push_back(centers[k]);
		polygon_mesh.color.push_back(color);
		for(int& corner: mycorners[k]) {
			polygon_mesh.position.push_back(corners[corner]);
			polygon_mesh.color.push_back(color);
		}

		// Generate triangle organization
		// The triangles are constructed as follows: a triangle is made up of
		// the center of a polygon, and the vertices of an edge of this polygon
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
	/// Initializes the shortest path algorithm
	/// assigning all corners INFINITY, except
	/// for the ocean polygons
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

	// Build the priority queue
	std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::function<bool(std::pair<int,int>, std::pair<int,int>)>> Q([] (const std::pair<int,int> &a, const std::pair<int,int> &b) {
   		return a.second > b.second;
	});
	for (int v = 0; v < N_corners; v++) {
		Q.push(std::pair<int,int>(v, corners[v].z));
	}
	
	while (!Q.empty()) {
		int u = Q.top().first;
		Q.pop();
		for (auto& adjacent: adjacents[u]) {
			int v = std::get<0>(adjacent);
			int poly_1 = std::get<1>(adjacent);
			int poly_2 = std::get<2>(adjacent);
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
}

/// Computes the smallest distance to a water source.
void scene_structure::compute_waterdists() {
	waterdists.resize(N);
	// Compute the source of the algorithm
	// and initializes the distances
	for (int idx = 0; idx < N; idx++) {
		if (biotopes[idx] == Biotope::Ocean || biotopes[idx] == Biotope::Lake) {
			waterdists[idx] = 0;
		} else {
			waterdists[idx] = INFINITY;
		}
	}

	// Build the priority queue
	std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::function<bool(std::pair<int,int>, std::pair<int,int>)>> Q([] (const std::pair<int,int> &a, const std::pair<int,int> &b) {
   		return a.second > b.second;
	});
	for (int v = 0; v < N; v++) {
		Q.push(std::pair<int,int>(v, waterdists[v]));
	}
	
	while (!Q.empty()) {
		int u = Q.top().first;
		Q.pop();
		for (auto& neighbor: neighbors[u]) {
			int v = std::get<0>(neighbor);
			float alt;
			if (biotopes[u] == Biotope::Ocean || biotopes[u] == Biotope::Lake) {
				alt = 0;
			} else {
				vec2 poly_u = { centers[u].x, centers[u].y };
				vec2 poly_v = { centers[v].x, centers[v].y };
				alt = waterdists[u] + norm(poly_u - poly_v);
			}

			if (alt < waterdists[v]) {
				waterdists[v] = alt;
				Q.push(std::pair<int,int>(v, waterdists[v]));
			}
		}
	}
}

/// Creates a new Voronoi structure with `n` clusters
/// and updates the class fields accordingly.
void scene_structure::create_voronoi(int n) {
	VoronoiDiagramGenerator vdg = VoronoiDiagramGenerator();
	std::vector<Point2>* sites = new std::vector<Point2>();
	BoundingBox bbox = BoundingBox(0, 10, 10, 0);

	// Create points, with possible duplicates
	// Use a temporary structure to hold data, we will remove duplicates afterwards
	// Heavily inspired by https://github.com/mdally/Voronoi/blob/master/examples/OpenGL_Example.cpp
	// in particular `genRandomSites()`
	std::vector<Point2> tmpSites;

	tmpSites.reserve(n);
	sites->reserve(n);

	Point2 s;

	for (unsigned int i = 0; i < n; ++i) {
		s.x = (rand() / (double)RAND_MAX)*10.0d;
		s.y = (rand() / (double)RAND_MAX)*10.0d;
		tmpSites.push_back(s);
	}

	// Remove any duplicates that exist
	// To do this, sort the temporary sites
	// and add them if they are not a duplicate.
	std::sort(tmpSites.begin(), tmpSites.end(), [] (const Point2& s1, const Point2& s2) {
		// Use lexicographic order to sort points
		if (s1.x < s2.x)
			return true;
		if (s1.x == s2.x && s1.y < s2.y)
			return true;

		return false;
	});
	sites->push_back(tmpSites[0]);
	for (Point2& s : tmpSites) {
		if (s != sites->back()) sites->push_back(s);
	}

	// We now know the number of Voronoi diagrams
	N = sites->size();
	centers.resize(N);
	biotopes.resize(N, Biotope::Land);
	neighbors.resize(N);
	mycorners.resize(N);

	Diagram* diagram = vdg.compute(*sites, bbox);

	// Store the centers
	int idx = 0;
	for (Cell* c : diagram->cells) {
		Point2& p = c->site.p;
		centers[idx] = {p.x, p.y, 0};
		idx++;
	}

	// Assign a number to each polygon and each vertex, and store it in
	// a map. These maps will help to find the number of each cell or vertex
	// in no time.
	std::map<Cell*, int> cells_ad;
	for (Cell* cell: diagram->cells) {
		int i = cells_ad.size();
		cells_ad.insert(std::pair<Cell*,int>(cell, i));
	}

	std::map<Point2*, int> vertices_ad;
	for (Point2* pt: diagram->vertices) {
		int i = corners.size();
		vec3 v = {pt->x, pt->y, 0.0f};
		corners.push_back(v);
		touches.push_back({});
		adjacents.push_back({});
		vertices_ad.insert(std::pair<Point2*,int>(pt, i));
	}

	N_corners = corners.size();

	// Store the corners and neighbors
	idx = 0;
	for (Cell* c: diagram->cells) {
		std::vector<std::tuple<int,int,int>> l;

        for (HalfEdge* halfedge: c->halfEdges) {
			Point2* vertA = halfedge->startPoint();
			Point2* vertB = halfedge->endPoint();
			vec3 v = {vertA->x, vertA->y, 0.0f};

			// For each point in the edge, retrieve its index in the corners list.
			// This assert and the following ones make sure that every vertex and
			// polygon must have been registered before in the first pass of the algo.
			auto found1 = vertices_ad.find(vertA);
			assert(found1 != vertices_ad.end());
			int pt1 = found1->second;
			// If there's some logic here, the list of corners is exactly
			// the list of the first vertex of each edge in the polygon
			mycorners[idx].push_back(pt1);

			v = {vertB->x, vertB->y, 0.0f};
			auto found2 = vertices_ad.find(vertB);
			assert(found2 != vertices_ad.end());
			int pt2 = found2->second;

			// Find the left polygon neighbor of the edge
			int neigh1;
			if (halfedge->edge->lSite) {
				auto found_idx1 = cells_ad.find(halfedge->edge->lSite->cell);
				assert(found_idx1 != cells_ad.end());
				neigh1 = found_idx1->second;
			} else {
				neigh1 = -1;
			}

			// Find the right polygon neighbor of the edge
			int neigh2;
			if (halfedge->edge->rSite) {
				auto found_idx2 = cells_ad.find(halfedge->edge->rSite->cell);
				assert(found_idx2 != cells_ad.end());
				neigh2 = found_idx2->second;
			} else {
				neigh2 = -1;
			}

			// Either neighbor 1 or neighbor2 of the edge is the current polygon
			assert(neigh1 == idx || neigh2 == idx && neigh1 != neigh2);

			// Which is why we can easily find the neighbor of this polygon
			// along this edge
			int idx_neigh = neigh1 == idx ? neigh2 : neigh1;

			// By convention, add to the adjacency list only if the index
			// of the neighboring polygon is smaller than the current one
			// (to avoid having the same adjdacents twice in the list). Takes
			// into account the case where the neighbor is -1 (no neighbor),
			// in which case we actually add to the adjacency list.
			if (idx_neigh < idx) {
				adjacents[pt1].push_back(std::make_tuple(pt2, neigh1, neigh2));
				adjacents[pt2].push_back(std::make_tuple(pt1, neigh1, neigh2));
			}

			// Since we are iterating over idx, we know
			// each idx will be added at most once.
			touches[pt1].push_back(idx);
			touches[pt2].push_back(idx);

			// Add the edge in the form of (neighbor, idx of 1st corner, idx of 2nd corner)
			l.push_back(std::make_tuple(idx_neigh, pt1, pt2));
		}
		neighbors[idx] = l;
		idx++;
	}

	// Always free up memory at the end
	delete sites;
	delete diagram;
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

/// Adds remaining biotopes.
void scene_structure::add_biotopes() {
	std::vector<float> orderedWaterdists = waterdists;
	sort(orderedWaterdists.begin(), orderedWaterdists.end());

	std::vector<float> orderedHeights = {};
	orderedHeights.resize(N);
	for (int i = 0; i < N; i++)
        orderedHeights[i] = corners[i].z;
	sort(orderedHeights.begin(), orderedHeights.end());

	// Find seperators between moisture and elevation zones
	float el2 = orderedHeights[N/4];
	float el3 = orderedHeights[N/2];
	float el4 = orderedHeights[3*N/4];
	float mois1 = orderedWaterdists[N/6];
	float mois2 = orderedWaterdists[2*N/6];
	float mois3 = orderedWaterdists[N/2];
	float mois4 = orderedWaterdists[2*N/3];
	float mois5 = orderedWaterdists[5*N/6];
	
	// We define snow biotopes based on elevation
    std::normal_distribution<double> distEl4(el4,SNOW_STD); 
    std::normal_distribution<double> distEl3(el3,SNOW_STD);
    std::normal_distribution<double> distEl2(el2,SNOW_STD);
	for (int idx = 0; idx < N; idx++) {
		if (biotopes[idx] != Biotope::Land)
			continue;
		if (centers[idx].z >= distEl4(randeng)) {
			if (waterdists[idx] <= mois3) biotopes[idx] = Biotope::Snow;
			else if (waterdists[idx] <= mois4) biotopes[idx] = Biotope::Tundra;
			else if (waterdists[idx] <= mois5) biotopes[idx] = Biotope::Bare;
			else biotopes[idx] = Biotope::Scorched;
		} else if (centers[idx].z >= distEl3(randeng)) {
			if (waterdists[idx] <= mois2) biotopes[idx] = Biotope::Taiga;
			else if (waterdists[idx] <= mois4) biotopes[idx] = Biotope::Shrubland;
			else biotopes[idx] = Biotope::TempDesert;
		} else if (centers[idx].z >= distEl2(randeng)) {
			if (waterdists[idx] <= mois1) biotopes[idx] = Biotope::TempRainForest;
			else if (waterdists[idx] <= mois3) biotopes[idx] = Biotope::TempDeciduousForest;
			else if (waterdists[idx] <= mois5) biotopes[idx] = Biotope::Grassland;
			else biotopes[idx] = Biotope::TempDesert;
		} else {
			if (waterdists[idx] <= mois2) biotopes[idx] = Biotope::TropicalRainForest;
			else if (waterdists[idx] <= mois4) biotopes[idx] = Biotope::TropicalSeasonalForest;
			else if (waterdists[idx] <= mois5) biotopes[idx] = Biotope::Grassland;
			else biotopes[idx] = Biotope::SubtropicalDesert;
		}
	}	
}