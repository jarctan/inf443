#include "scene.hpp"

const float MOIS_STD = 0.2f;
const float EL_STD = 0.2f;
/// Number of relaxations to perform on the Voronoi diagram.
const int RELAX_CNT = 2;
/// (Estimate) number of polygons to use to generate the diagram.
const int POLYGON_CNT = 30000;
/// SIze of the generated terrain (height and width)
const float SIZE = 20.0f;

mesh cylinder(float r, float h);

void scene_structure::handleKeyPress(GLFWwindow* window, int key, int action) {

	inputs_keyboard_parameters const& keyboard = inputs.keyboard;
	
	//block handling activation and deactivation of camera movement
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		cameraCanMove = !cameraCanMove;
		if (cameraCanMove) {
			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
		}
		else {
			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
		}
	}

	//block handling camera movement
	if (cameraCanMove) {
		camera_head& camera = environment.camera;

		float const dt = timer.update();
		vec3 displacement;

		if (key == GLFW_KEY_W && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			vec3 frontDirection = vec3(camera.front().x, camera.front().y, 0);
			frontDirection = frontDirection / norm(frontDirection);
			displacement = speed * dt * frontDirection;
		}
		else if (key == GLFW_KEY_A && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			vec3 leftDirection = vec3(-camera.front().y, camera.front().x, 0);
			leftDirection = leftDirection / norm(leftDirection);
			displacement = speed * dt * leftDirection;
		}
		else if (key == GLFW_KEY_D && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			vec3 rightDirection = vec3(camera.front().y, -camera.front().x, 0);
			rightDirection = rightDirection / norm(rightDirection);
			displacement = speed * dt * rightDirection;
		}
		else if (key == GLFW_KEY_S && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			vec3 backDirection = vec3(-camera.front().x, -camera.front().y, 0);
			backDirection = backDirection / norm(backDirection);
			displacement = speed * dt * backDirection;
		}
		else if (key == GLFW_KEY_SPACE && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			displacement = speed * dt * vec3(0, 0, 1);
		}
		else if (key == GLFW_KEY_ENTER && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			displacement = speed * dt * vec3(0, 0, -1);
		}

		camera.position_camera += displacement;
	}
}

void scene_structure::handleMouseMove(GLFWwindow* window) {
	if (cameraCanMove) {
		int screenWidth = 0, screenHeight = 0;
		glfwGetWindowSize(window, &screenWidth, &screenHeight);
		double xpos, ypos;
		glfwGetCursorPos(window, &xpos, &ypos);
		glfwSetCursorPos(window, screenWidth / 2.0f, screenHeight / 2.0f);
		camera_head& camera = environment.camera;
		float const dt = timer.update();

		rotation_transform initialOrientation = rotation_transform::from_axis_angle(vec3(1, 0, 0), initial_camera_pitch);

		camera_yaw += (screenWidth / 2.0f - xpos) * dt * mouseSpeed;
		camera_pitch += (screenHeight / 2.0f - ypos) * dt * mouseSpeed;

		rotation_transform r_pitch = rotation_transform::from_axis_angle(vec3(1, 0, 0), camera_pitch);
		rotation_transform r_yaw = rotation_transform::from_axis_angle(vec3(0, 1, 0), camera_yaw);

		camera.orientation_camera = initialOrientation * r_yaw * r_pitch;
	}
}

/// This function is called only once at the beginning of the program
/// and initializes the meshes, diagrams and other structures.
void scene_structure::initialize() {
	// Set the behavior of the camera and its initial position
	randeng = default_random_engine(time(NULL));

	auto start = high_resolution_clock::now();
	create_voronoi(POLYGON_CNT);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << "Voronoi diagram in " << duration.count() << "ms [OK]" << endl;

	start = high_resolution_clock::now();
	create_ocean_border();
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	cout << "Ocean border in " << duration.count() << "ms [OK]" << endl;

	start = high_resolution_clock::now();
	compute_heights();
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	cout << "Compute height in " << duration.count() << "ms [OK]" << endl;

	start = high_resolution_clock::now();
	compute_waterdists();
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	cout << "Water distance in " << duration.count() << "ms [OK]" << endl;

	// Compute the elevation of the centers of the polygons
	start = high_resolution_clock::now();
	smooth_centers();
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	cout << "Smooth centers in " << duration.count() << "ms [OK]" << endl;

	start = high_resolution_clock::now();
	add_biotopes();
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	cout << "Biotopes in " << duration.count() << "ms [OK]" << endl;

	start = high_resolution_clock::now();
	add_wind();
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	cout << "Wind field in " << duration.count() << "ms [OK]" << endl;

	start = high_resolution_clock::now();
	create_terrain();
	stop = high_resolution_clock::now();
	cout << "Terrain in " << duration.count() << "ms [OK]" << endl;

	// Initialize the skybox
	// The path must contain 6 texture images
	skybox.initialize("assets/skybox/");

    float r = 0.1f;
    float h = 0.3f;
	mesh cylinder_mesh = cylinder(r, h);
	ship.initialize(cylinder_mesh, "Ship");
	for (int idx = 0; idx < N; idx++) {
		if (biotopes[idx] == Biotope::Ocean) {
			ship.transform.translation = centers[idx];
			break;
		}
	}

	timer.scale = 0.5f;

	// Create a visual frame representing the coordinate system
	global_frame.initialize(mesh_primitive_frame(), "Frame");
	environment.camera.position_camera = { 5.0f, 5.0f, 10.0f };
	environment.camera.manipulator_rotate_roll_pitch_yaw(0, camera_pitch, camera_yaw); //initial rotation value

	timer.scale = 0.5f;
}

// This function is constantly called at every frame
void scene_structure::display() {
	// Display of the skybox
	// Note: The skybox must be drawn before the other shapes 
	// Skybox is displayed without writing in the z-buffer.
	// In displaying it first, the cube appears beyond any other shape.
	draw(skybox, environment);

	// Set the light to the current position of the camera
	environment.light = environment.camera.position();

	// conditional display of the global frame (set via the GUI)
	if (gui.display_frame)
		draw(global_frame, environment);

	// Update the current time
	timer.update();
	auto wind = windfield[(int) ship.transform.translation.x][(int) ship.transform.translation.y];
	float ship_x = ship.transform.translation.x + wind.first;
	float ship_y = ship.transform.translation.y + wind.second;
	if (ship_x >= 0.0f && ship_y >= 0.0f && ship_x <= (float) SIZE && ship_y <= (float) SIZE) {
		ship.transform.translation.x = ship_x;
		ship.transform.translation.y = ship_y;
		ship.transform.translation.z = heightfield[(int) ship_x][(int) ship_y];
	}

	// Draw the lines between the clusters, and the borders
	draw(terrain, environment);
	draw(ship, environment);
}

/// Displays the GUI elements.
void scene_structure::display_gui() {
	ImGui::Checkbox("Frame", &gui.display_frame);
}

/// Creates a terrain mesh based on the Voronoi diagram.
void scene_structure::create_terrain() {
	// Initializes heightmap
	heightfield.resize((int) SIZE + 1, {});
	for (int i = 0; i < (int) SIZE + 1; i++) {
		heightfield[i].resize((int) SIZE + 1, INFINITY);
	}

    // Terrain geometry
	mesh terrain_mesh;
    for (int k = 0; k < N; k++) {
		// The following convention is being used:
		// the first point is the centers of the polygon,
		// then the following points are its corners.
		mesh polygon_mesh;
		const int center_idx = 0;

		vec3 color;
		if (biotopes[k] == Biotope::Ocean)
			color = vec3(0.0f,0.412f,0.58f);
		else if (biotopes[k] == Biotope::Lake)
			color = vec3(0.2f,0.59f,0.885f);
		else if (biotopes[k] == Biotope::Snow)
			color = vec3(1.0f,1.0f,1.0f);
		else if (biotopes[k] == Biotope::Tundra)
			color = vec3(0.867f,0.867f,0.733f);
		else if (biotopes[k] == Biotope::Bare)
			color = vec3(0.733f,0.733f,0.733f);
		else if (biotopes[k] == Biotope::Scorched)
			color = vec3(0.6f,0.6f,0.6f);
		else if (biotopes[k] == Biotope::Taiga)
			color = vec3(0.8f,0.831f,0.733f);
		else if (biotopes[k] == Biotope::Shrubland)
			color = vec3(0.769f,0.8f,0.733f);
		else if (biotopes[k] == Biotope::TempDesert)
			color = vec3(0.894f,0.91f,0.792f);
		else if (biotopes[k] == Biotope::TempRainForest)
			color = vec3(0.643f,0.769f,0.659f);
		else if (biotopes[k] == Biotope::TempDeciduousForest)
			color = vec3(0.706f,0.788f,0.663f);
		else if (biotopes[k] == Biotope::Grassland)
			color = vec3(0.769f,0.831f,0.667f);
		else if (biotopes[k] == Biotope::TropicalRainForest)
			color = vec3(0.612f,0.733f,0.663f);
		else if (biotopes[k] == Biotope::TropicalSeasonalForest)
			color = vec3(0.663f,0.8f,0.643f);
		else if (biotopes[k] == Biotope::SubtropicalDesert)
			color = vec3(0.914f,0.867f,0.78f);
		else
			color = vec3(0.6f,0.85f,0.5f);

		heightfield[(int) centers[k].x][(int) centers[k].y] = centers[k].z;
		polygon_mesh.position.push_back(centers[k]);
		polygon_mesh.color.push_back(color);
		for(int& corner: mycorners[k]) {
			heightfield[(int) corners[corner].x][(int) corners[corner].y] = corners[corner].z;
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
	priority_queue<pair<int,int>, vector<pair<int,int>>, function<bool(pair<int,int>, pair<int,int>)>> Q([] (const pair<int,int> &a, const pair<int,int> &b) {
   		return a.second > b.second;
	});
	for (int v = 0; v < N_corners; v++) {
		Q.push(pair<int,int>(v, corners[v].z));
	}
	
	while (!Q.empty()) {
		int u = Q.top().first;
		Q.pop();
		for (auto& adjacent: adjacents[u]) {
			int v = adjacent.vertex;
			int poly_1 = adjacent.polyA;
			int poly_2 = adjacent.polyB;
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
				vec2 pt = normalize((corner_u + corner_v) / 2);
				float noise = noise_perlin(pt, 2, 0.291f, 10.0f)/15.0f;
				// float noise = std::pow(noise_perlin(pt, 6, 0.291f, 5.268f),1.7f); when multiplied by norm
				// float noise = noise_perlin(pt, 6, 0.507f, 1.982f); when multiplied by norm
				dist = noise;
			}

			float alt = max(corners[u].z + dist, 0.0f);
			if (alt < corners[v].z) {
				corners[v].z = alt;
				Q.push(pair<int,int>(v, corners[v].z));
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
	priority_queue<pair<int,int>, vector<pair<int,int>>, function<bool(pair<int,int>, pair<int,int>)>> Q([] (const pair<int,int> &a, const pair<int,int> &b) {
   		return a.second > b.second;
	});
	for (int v = 0; v < N; v++) {
		Q.push(pair<int,int>(v, waterdists[v]));
	}
	
	while (!Q.empty()) {
		int u = Q.top().first;
		Q.pop();
		for (auto& neighbor: neighbors[u]) {
			int v = neighbor.polygon;
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
				Q.push(pair<int,int>(v, waterdists[v]));
			}
		}
	}
}

/// Creates a new Voronoi structure with `n` clusters
/// and updates the class fields accordingly.
void scene_structure::create_voronoi(int n) {
	VoronoiDiagramGenerator vdg = VoronoiDiagramGenerator();
	vector<Point2>* sites = new vector<Point2>();
	BoundingBox bbox = BoundingBox(0, SIZE, SIZE, 0);

	// Create points, with possible duplicates
	// Use a temporary structure to hold data, we will remove duplicates afterwards
	// Heavily inspired by https://github.com/mdally/Voronoi/blob/master/examples/OpenGL_Example.cpp
	// in particular `genRandomSites()`
	vector<Point2> tmpSites;

	tmpSites.reserve(n);
	sites->reserve(n);

	Point2 s;

	for (unsigned int i = 0; i < n; ++i) {
		s.x = (rand() / (double)RAND_MAX)*(double) SIZE;
		s.y = (rand() / (double)RAND_MAX)*(double) SIZE;
		tmpSites.push_back(s);
	}

	// Remove any duplicates that exist
	// To do this, sort the temporary sites
	// and add them if they are not a duplicate.
	sort(tmpSites.begin(), tmpSites.end(), [] (const Point2& s1, const Point2& s2) {
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

	// Add relaxation (place gems at the centroid
	// of their polygons and recompute Voronoi diagrams)
	for (int i = 0; i < RELAX_CNT; i++) {
		Diagram* prevDiagram = diagram;
		diagram = vdg.relax();
		delete prevDiagram;
	}

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
	map<Cell*, int> cells_ad;
	for (Cell* cell: diagram->cells) {
		int i = cells_ad.size();
		cells_ad.insert(pair<Cell*,int>(cell, i));
	}

	map<Point2*, int> vertices_ad;
	for (Point2* pt: diagram->vertices) {
		int i = corners.size();
		vec3 v = {pt->x, pt->y, 0.0f};
		corners.push_back(v);
		touches.push_back({});
		adjacents.push_back({});
		vertices_ad.insert(pair<Point2*,int>(pt, i));
	}

	N_corners = corners.size();

	// Store the corners and neighbors
	idx = 0;
	for (Cell* c: diagram->cells) {
		vector<Neighbor> l;

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
				neigh1 = 0;
			}

			// Find the right polygon neighbor of the edge
			int neigh2;
			if (halfedge->edge->rSite) {
				auto found_idx2 = cells_ad.find(halfedge->edge->rSite->cell);
				assert(found_idx2 != cells_ad.end());
				neigh2 = found_idx2->second;
			} else {
				neigh2 = 0;
			}

			// Either neighbor 1 or neighbor2 of the edge is the current polygon
			assert((neigh1 == idx || neigh2 == idx) && neigh1 != neigh2);

			// Which is why we can easily find the neighbor of this polygon
			// along this edge
			int idx_neigh = neigh1 == idx ? neigh2 : neigh1;

			// By convention, add to the adjacency list only if the index
			// of the neighboring polygon is smaller than the current one
			// (to avoid having the same adjdacents twice in the list). Takes
			// into account the case where the neighbor is -1 (no neighbor),
			// in which case we actually add to the adjacency list.
			if (idx_neigh < idx) {
				adjacents[pt1].push_back({pt2, neigh1, neigh2});
				adjacents[pt2].push_back({pt1, neigh1, neigh2});
			}

			// Since we are iterating over idx, we know
			// each idx will be added at most once.
			touches[pt1].push_back(idx);
			touches[pt2].push_back(idx);

			// Add the edge in the adjacency list
			l.push_back({idx_neigh, pt1, pt2});
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
	// Random seed computed once to generate different perlin noise
	// at each time. This seed is used as an offset of the position
	// to determine where to look in the (bidimentional, infinite) texture.
	float delta_x = (rand() / (double)RAND_MAX)*20.0f;
	float delta_y = (rand() / (double)RAND_MAX)*20.0f;

	for(int idx = 0; idx < N; idx++) {
		for (auto& corner_idx: mycorners[idx]) {
			vec3& corner = corners[corner_idx];
			if (corner.x <= 0.2f || corner.y <= 0.2f || corner.x >= (float) SIZE - 0.2f || corner.y  >= (float) SIZE - 0.2f) {
				  biotopes[idx] = Biotope::Ocean;
			}
		}
		vec2 pt = {centers[idx].x / (float) 10.0f + delta_x, centers[idx].y / (float) 10.0f + delta_y};
		// First proposition: if (noise_perlin(pt, 6, 0.291f, 5.268f) > 0.6f) {
		if (noise_perlin(pt, 6, 0.291f, 5.268f) > 0.8f // Ocean
			|| noise_perlin({pt.x * 1.5f, pt.y * 1.5f}, 2, 0.355f, 3.603f) > 1.1f // Lakes
		) {
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
			mean += corners[adjacent.vertex].z;
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
	vector<float> orderedWaterdists = waterdists;
	sort(orderedWaterdists.begin(), orderedWaterdists.end());

	vector<float> orderedHeights = {};
	orderedHeights.resize(N);
	for (int i = 0; i < N; i++)
        orderedHeights[i] = corners[i].z;
	sort(orderedHeights.begin(), orderedHeights.end());

	// Find highest, finite height
	float max_height;
	for (int i = N-1; i >= 0; i--) {
		if (orderedHeights[i] < INFINITY) {
			max_height = orderedHeights[i];
			break;
		}
	}

	// Find seperator values between elevation zones
	float el2 = max_height/4;
	float el3 = max_height/2;
	float el4 = 3*max_height/4;

	// Find separator values between moisture zones
	float mois1 = orderedWaterdists[N/6];
	float mois2 = orderedWaterdists[2*N/6];
	float mois3 = orderedWaterdists[N/2];
	float mois4 = orderedWaterdists[2*N/3];
	float mois5 = orderedWaterdists[5*N/6];
	
	// We define snow biotopes based on elevation
    normal_distribution<double> distEl4(el4,EL_STD); 
    normal_distribution<double> distEl3(el3,EL_STD);
    normal_distribution<double> distEl2(el2,EL_STD);
    normal_distribution<double> distMois1(mois1,MOIS_STD); 
    normal_distribution<double> distMois2(mois2,MOIS_STD);
    normal_distribution<double> distMois3(mois3,MOIS_STD);
    normal_distribution<double> distMois4(mois4,MOIS_STD);
    normal_distribution<double> distMois5(mois5,MOIS_STD);
	for (int idx = 0; idx < N; idx++) {
		if (biotopes[idx] != Biotope::Land)
			continue;
		if (centers[idx].z >= distEl4(randeng)) {
			if (waterdists[idx] <= distMois3(randeng)) biotopes[idx] = Biotope::Snow;
			else if (waterdists[idx] <= distMois4(randeng)) biotopes[idx] = Biotope::Tundra;
			else if (waterdists[idx] <= distMois5(randeng)) biotopes[idx] = Biotope::Bare;
			else biotopes[idx] = Biotope::Scorched;
		} else if (centers[idx].z >= distEl3(randeng)) {
			if (waterdists[idx] <= distMois2(randeng)) biotopes[idx] = Biotope::Taiga;
			else if (waterdists[idx] <= distMois4(randeng)) biotopes[idx] = Biotope::Shrubland;
			else biotopes[idx] = Biotope::TempDesert;
		} else if (centers[idx].z >= distEl2(randeng)) {
			if (waterdists[idx] <= distMois1(randeng)) biotopes[idx] = Biotope::TempRainForest;
			else if (waterdists[idx] <= distMois3(randeng)) biotopes[idx] = Biotope::TempDeciduousForest;
			else if (waterdists[idx] <= distMois5(randeng)) biotopes[idx] = Biotope::Grassland;
			else biotopes[idx] = Biotope::TempDesert;
		} else {
			if (waterdists[idx] <= distMois2(randeng)) biotopes[idx] = Biotope::TropicalRainForest;
			else if (waterdists[idx] <= distMois4(randeng)) biotopes[idx] = Biotope::TropicalSeasonalForest;
			else if (waterdists[idx] <= distMois5(randeng)) biotopes[idx] = Biotope::Grassland;
			else biotopes[idx] = Biotope::SubtropicalDesert;
		}
	}	
}


/// Adds wind.
void scene_structure::add_wind() {
	for (int i = 0; i < (int) SIZE; i++) {
		windfield.push_back({});
		for (int j = 0; j < (int) SIZE; j++) {
			windfield[i].push_back(pair<float,float>(0.005f,-0.005f));
		}
	}
}

/// Create a cylinder
/// TODO: change or remove to create a real ship.
mesh cylinder(float r, float h) {
    // Number of samples of the terrain is N x N
    int N = 20;

    mesh cylinder; // temporary terrain storage (CPU only)
    cylinder.position.resize(N*N);
    cylinder.uv.resize(N*N);

    // Fill terrain geometry
    for(int ku=0; ku<N; ++ku)
    {
        for(int kv=0; kv<N; ++kv)
        {
            // Compute local parametric coordinates (u,v) \in [0,1]
            float u = ku/(N-1.0f);
            float v = kv/(N-1.0f);

            // Compute the local surface function
            vec3 p = {r*cos(2* Pi *u), r*sin(2* Pi *u), h*(v-0.5f)};

            // Store vertex coordinates
            cylinder.position[kv+N*ku] = p;
            cylinder.uv[kv+N*ku] = {u, v};
        }
    }

    // Generate triangle organization
    for(int ku=0; ku<N-1; ++ku)
    {
        for(int kv=0; kv<N-1; ++kv)
        {
            int idx = kv + N*ku;

            uint3 triangle_1 = {idx, idx+1+N, idx+1};
            uint3 triangle_2 = {idx, idx+N, idx+1+N};

            cylinder.connectivity.push_back(triangle_1);
            cylinder.connectivity.push_back(triangle_2);
        }
    }

	cylinder.fill_empty_field();
    return cylinder;
}