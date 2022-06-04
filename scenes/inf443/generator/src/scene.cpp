#include "scene.hpp"

scene_structure::~scene_structure() {
	for (Windsource* source: windsources) {
		delete source;
	}
}

void scene_structure::handleKeyPress(GLFWwindow* window, int key, int action) {
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
	cout << "Terrain generated in " << duration.count() << "ms [OK]" << endl;

	start = high_resolution_clock::now();
	add_ships();
	stop = high_resolution_clock::now();
	cout << "Ships added in " << duration.count() << "ms [OK]" << endl;

	start = high_resolution_clock::now();
	add_snowflakes();
	stop = high_resolution_clock::now();
	cout << "Snow added in " << duration.count() << "ms [OK]" << endl;

	start = high_resolution_clock::now();
	add_birds();
	stop = high_resolution_clock::now();
	cout << "Birds added in " << duration.count() << "ms [OK]" << endl;

	start = high_resolution_clock::now();
	add_cloud();
	stop = high_resolution_clock::now();
	cout << "Cloud added in " << duration.count() << "ms [OK]" << endl;

	// Initialize the skybox
	// The path must contain 6 texture images
	skybox.initialize("assets/skybox/");

	// Create a visual frame representing the coordinate system
	global_frame.initialize(mesh_primitive_frame(), "Frame");
	environment.camera.position_camera = { 5.0f, 5.0f, 10.0f };
	environment.camera.manipulator_rotate_roll_pitch_yaw(0, camera_pitch, camera_yaw); //initial rotation value

	timer.scale = 0.5;
	camera_timer.scale = 0.1;
	movement_timer.scale = 0.1;
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
	float dt = timer.update();
	float t = timer.t;

	// Update the wind primitives
	for (Windsource* source : windsources) {
		source->step(timer.scale);
	}

	// Draw the terrain and the sea
	draw(terrain, environment);
	draw(sea, environment);
	
	// Draw birds
	update_birds_positions(dt, n_birds, birds);
	for (int i = 0; i < n_birds; i++) {
		birds[i].bird_drawable["body"].transform.translation = vec3( 0,0.2f * std::cos(5 * t),0 ) + birds[i].pos;
		birds[i].bird_drawable["wing_L"].transform.rotation = rotation_transform::from_axis_angle({ 0,0,1 }, Pi / 8.0f * (1 + std::cos(20 * t)));
		birds[i].bird_drawable["wing_R"].transform.rotation = rotation_transform::from_axis_angle({ 0,0,-1 }, Pi / 8.0f * (1 + std::cos(20 * t)));
		birds[i].bird_drawable["head"].transform.translation = vec3{ 0,0,1.8 } + vec3{ 0,0,0.1f * std::cos(2.5 * t) };
		birds[i].bird_drawable["head"].transform.rotation = rotation_transform::from_axis_angle({ 0,1,0 }, Pi / 8.0f * std::cos(5 * t));
		
		birds[i].bird_drawable.update_local_to_global_coordinates();

		birds[i].bird_drawable["body"].transform.rotation = rotation_transform::between_vector({0.0f, 1.0f, 0.0f}, normalize(birds[i].speed));
		draw(birds[i].bird_drawable, environment);
	}

	// Compute the wind direction for each ship
	// Taking into account the influence of each ship on
	// the other ones
	vector<vec2> wind_dir;
	wind_dir.resize(ships.size());
	for (int i = 0; i < ships.size(); i++) {
		Ship& ship = ships[i];
		if (ship.outofscope) continue;

		vec2 wind;
		for (Windsource* source: windsources) {
			wind += source->field(ship.pos);
		}
		for (Ship& other_ship: ships) {
			Revultion source = Revultion(other_ship.pos, 0.3f, 0.3f, 0.002f);
			wind += source.field(ship.pos);
		}

		wind_dir[i] = wind;
		if (norm(wind) > 0.0001f) {
			ship.dir = normalize(wind);
		}
	}

	// Update the ships according to the previously
	// computed wind direction
	// Mark if necessary newly out of scope ships
	for (int i = 0; i < ships.size(); i++) {
		Ship& ship = ships[i];
		if (ship.outofscope) continue;

		vec2 new_pos = ship.pos + wind_dir[i];

		// If the ship is still in the terrain, update its position and draw it
		if (new_pos.x >= 0.0f && new_pos.y >= 0.0f && new_pos.x <= (float) TERRAIN_SIZE && new_pos.y <= (float) TERRAIN_SIZE) {
			// Update if necessary the polygon where the ships stands
			// Since we rely on Voronoi cells, it is known that if you are closer
			// to the center of a neighboring centroid than to the center of your
			// own polygon, you now belong to the neighboring polygon
			int new_polygon = ship.polygon;
			vec2 my_center = { centers[ship.polygon].x, centers[ship.polygon].y };
			for(Neighbor& neighbor: neighbors[ship.polygon]) {
				vec2 new_center = { centers[neighbor.polygon].x, centers[neighbor.polygon].y };
				if (norm(ship.pos-new_center) < norm(ship.pos-my_center)) {
					new_polygon = neighbor.polygon;
					break;
				}
			}


			// Accept the updates only if we do NOT land on the shores
			if (biotopes[new_polygon] == Biotope::Ocean) {
				ship.pos = new_pos;
				ship.polygon = new_polygon;
			}

			ship_drawable["ship"].transform.translation = { ship.pos.x, ship.pos.y, 0 };
			// The first reorientation is needed because the OBJ is not correcly oriented along this axis
			ship_drawable["ship"].transform.rotation = rotation_transform::from_axis_angle({ 1,0,0 }, M_PI/2);

			// We orient the ship
			ship_drawable["ship"].transform.rotation = rotation_transform::between_vector({ 0.0f, -1.0f, 0.0f }, {ship.dir.x, ship.dir.y, 0.0f }) * ship_drawable["ship"].transform.rotation;
			
			// We place the sail at the same position as the ship
			ship_drawable["sail"].transform = ship_drawable["ship"].transform;

			// We update the coordinates of the whole hierarchy
			ship_drawable.update_local_to_global_coordinates();
	
			draw(ship_drawable, environment);
		} else {
			ship.outofscope = true;
		}
	}

	remove_if(ships.begin(), ships.end(), [] (Ship& ship) {
		return ship.outofscope;
	});

	display_semi_transparent();
}

void scene_structure::display_semi_transparent() {
	// Enable use of alpha component as color blending for transparent elements
	//  alpha = current_color.alpha
	//  new color = previous_color * alpha + current_color * (1-alpha)
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Disable depth buffer writing
	//  - Transparent elements cannot use depth buffer
	//  - They are supposed to be display from furest to nearest elements
	glDepthMask(false);

	// For each particle of the cloud, update its position
	// according to the wind direction and draw it
	for (vec3& particle: particles) {
		vec2 wind;
		vec2 pos2D = { particle.x, particle.y };
		for (Windsource* source: windsources) {
			wind += source->field(pos2D);
		}
		pos2D += wind;
		particle.x = pos2D.x;
		particle.y = pos2D.y;

		// Use these guiding particles to draw many particles around
		// Seems however to be really resources-intensive, so we may have to limit the number of particles drawn
		// to only one. If you computer can handle it, you can increase PARTICLES_PER_GUIDE and see the result.
		// We use a slight perlin noise to have a deterministic, continuous noise
		for (int i = -PARTICLES_PER_GUIDE; i <= PARTICLES_PER_GUIDE; i++) {
			for (int j = -PARTICLES_PER_GUIDE; j <= PARTICLES_PER_GUIDE; j++) {
				float perlin_x = noise_perlin({particle.x / (float) TERRAIN_SIZE, particle.y / (float) TERRAIN_SIZE}, 4, 3.26f, 2.268f) / 500;
				float perlin_y = noise_perlin({particle.x / (float) TERRAIN_SIZE, particle.y / (float) TERRAIN_SIZE}, 4, 3.26f, 2.978f) / 500;
				cloud.transform.translation = { particle.x + i * PARTICLE_SIZE + perlin_x, particle.y + perlin_y + j * PARTICLE_SIZE, particle.z } ;
				draw(cloud, environment);
			}
		}
	}

	// For each snowflake, update and print
	for (Snowflake& particle: snowflakes) {
		vec2 wind;
		vec2 pos2D = { particle.pos.x, particle.pos.y };
		for (Windsource* source: windsources) {
			wind += source->field(pos2D);
		}
		pos2D += wind;
		particle.pos.x = pos2D.x;
		particle.pos.y = pos2D.y;
		particle.speed_z -= 0.0000981;
		particle.pos.z += particle.speed_z * timer.scale;

		// Find
		// Update if necessary the polygon where the snowflake stands	
		// Since we rely on Voronoi cells, it is known that if you are closer
		// to the center of a neighboring centroid than to the center of your
		// own polygon, you now belong to the neighboring polygon
		vec2 my_center = { centers[particle.polygon].x, centers[particle.polygon].y };
		vec2 my_pos = { particle.pos.x, particle.pos.y };
		for(Neighbor& neighbor: neighbors[particle.polygon]) {
			vec2 new_center = { centers[neighbor.polygon].x, centers[neighbor.polygon].y };
			if (norm(my_pos-new_center) < norm(my_pos-my_center)) {
				particle.polygon = neighbor.polygon;
				break;
			}
		}

		// If we reach the ground, spawn again in the sky
		if (particle.pos.z < centers[particle.polygon].z) {
			particle.pos.x = centers[particle.initial_polygon].x;
			particle.pos.y = centers[particle.initial_polygon].y;
			particle.pos.z = snowheight(randeng);
			particle.speed_z = 0;

			// If there snow on a non-snow biotope, change the biotope
			// on a random basis
			if (biotopes[particle.polygon] != Biotope::Snow 
				&& biotopes[particle.polygon] != Biotope::Lake
				&& biotopes[particle.polygon] != Biotope::Ocean
				&& rand() % 75 == 0) {
				cout << "New snow biotope added" << endl;
				biotopes[particle.polygon] = Biotope::Snow;
				terrain.clear();
				sea.clear();
				create_terrain();
			}
		}

		snowflake.transform.translation = particle.pos;
		draw(snowflake, environment);
	}

	// Don't forget to re-activate the depth-buffer write
	glDepthMask(true);
	glDisable(GL_BLEND);
}

/// Displays the GUI elements.
void scene_structure::display_gui() {
	ImGui::Checkbox("Frame", &gui.display_frame);
}

/// Creates a terrain mesh based on the Voronoi diagram.
void scene_structure::create_terrain() {
    // Terrain geometry
	mesh terrain_mesh;
    // Sea geometry
	mesh sea_mesh;
    for (int k = 0; k < N; k++) {
		// The following convention is being used:
		// the first point is the centers of the polygon,
		// then the following points are its corners.
		mesh polygon_mesh;
		const int center_idx = 0;

		vec3 color;
		if (biotopes[k] == Biotope::Ocean || biotopes[k] == Biotope::Shore)
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

		if (biotopes[k] == Biotope::Ocean || biotopes[k] == Biotope::Shore)
			sea_mesh.push_back(polygon_mesh);
		else
			terrain_mesh.push_back(polygon_mesh);
	}

	// Initialize and sets the right color for the land
	terrain.initialize(terrain_mesh, "Terrain");
	terrain.shading.phong.diffuse = 0.2f;
	terrain.shading.phong.ambient = 0.8f;
	terrain.shading.phong.specular = 0.0f; // non-specular land material

	// Initialize and sets the right color for the land
	sea.initialize(sea_mesh, "Sea");
	sea.shading.phong.diffuse = 0.0f;
	sea.shading.phong.ambient = 0.8f;
	sea.shading.phong.specular = 0.0f; // non-specular land material
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
			if (poly_1 != -1 && biotopes[poly_1] == Biotope::Ocean) {
				if (poly_2 != -1 && biotopes[poly_2] == Biotope::Lake) biotopes[poly_2] = Biotope::Ocean;
				dist = 0;
			} else if (poly_2 != -1 && biotopes[poly_2] == Biotope::Ocean) {
				if (poly_1 != -1 && biotopes[poly_1] == Biotope::Lake) biotopes[poly_1] = Biotope::Ocean;
				dist = 0;
			} else if ((poly_1 != -1 && biotopes[poly_1] == Biotope::Lake) || (poly_2 != -1 && biotopes[poly_2] == Biotope::Lake)) {
				dist = 0;
			} else {
				vec2 corner_u = { corners[u].x, corners[u].y };
				vec2 corner_v = { corners[v].x, corners[v].y };
				vec2 pt = normalize((corner_u + corner_v) / 2);
				float noise = noise_perlin(pt, 2, 0.291f, 10.0f)/15.0f;
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
			if (v == -1) continue;
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
	BoundingBox bbox = BoundingBox(0, TERRAIN_SIZE, TERRAIN_SIZE, 0);

	// Create points, with possible duplicates
	// Use a temporary structure to hold data, we will remove duplicates afterwards
	// Heavily inspired by https://github.com/mdally/Voronoi/blob/master/examples/OpenGL_Example.cpp
	// in particular `genRandomSites()`
	vector<Point2> tmpSites;

	tmpSites.reserve(n);
	sites->reserve(n);

	Point2 s;

	for (unsigned int i = 0; i < n; ++i) {
		s.x = (rand() / (double)RAND_MAX)*(double) TERRAIN_SIZE;
		s.y = (rand() / (double)RAND_MAX)*(double) TERRAIN_SIZE;
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
	float delta_x = (rand() / (double)RAND_MAX)*1000.0f;
	float delta_y = (rand() / (double)RAND_MAX)*1000.0f;

	for(int idx = 0; idx < N; idx++) {
		for (auto& corner_idx: mycorners[idx]) {
			vec3& corner = corners[corner_idx];
			if (corner.x <= 0.2f || corner.y <= 0.2f || corner.x >= (float) TERRAIN_SIZE - 0.2f || corner.y  >= (float) TERRAIN_SIZE - 0.2f) {
				  biotopes[idx] = Biotope::Ocean;
			}
		}
		vec2 pt = {centers[idx].x / (float) 10.0f + delta_x, centers[idx].y / (float) 10.0f + delta_y};
		// First proposition: if (noise_perlin(pt, 6, 0.291f, 5.268f) > 0.6f) {
		if (noise_perlin(pt, 6, 0.291f, 5.268f) > 0.8f // Ocean
			|| noise_perlin({pt.x * 1.5f, pt.y * 1.5f}, 2, 0.355f, 3.603f) > 1.1f // Lakes
			|| norm(vec2 { centers[idx].x, centers[idx].y } - vec2 {(float) TERRAIN_SIZE / 2, (float) TERRAIN_SIZE / 2}) >= (0.75f + sin(7*atan(centers[idx].y/centers[idx].x))*0.15f) * (float) TERRAIN_SIZE / 1.4142f
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
			biotopes[idx] = Biotope::Snow;
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
	// Find the shore biotopes (these are
	// the sea biotopes that are near the land)
	for (int idx = 0; idx < N; idx++) {
		if (biotopes[idx] == Biotope::Ocean) {
			for (Neighbor& neighbor: neighbors[idx]) {
				int poly = neighbor.polygon;
				if (poly != -1 && biotopes[poly] != Biotope::Ocean) {
					biotopes[idx] = Biotope::Shore;
				}
			}
		}
	}

	// Add cyclones
	for (int i = 0; i < 5; i++) {
		Cyclone* cyclone = new Cyclone(vec2 { (rand() / (float) RAND_MAX)*(float) TERRAIN_SIZE, (rand() / (float) RAND_MAX)*(float) TERRAIN_SIZE}, 0.002f, 1.5f, { 0.001, 0.001 }, rand() % 2 == 0);
		windsources.push_back(cyclone);
	}

	// Add gusts
	for (int i = 0; i < 5; i++) {
		Gust* gust = new Gust(vec2 { (rand() / (float) RAND_MAX)*(float) TERRAIN_SIZE, (rand() / (float) RAND_MAX)*(float) TERRAIN_SIZE}, 0.001f * normalize(vec2 { rand(), rand() }), 0.002f, 2);
		windsources.push_back(gust);
	}

	// And a final one for the clouds
	Gust* gust = new Gust({0.0f, 0.0f}, {0.001f, 0.001f}, 0.005f, 1.0f);
	windsources.push_back(gust);
}

/// Adds ships.
void scene_structure::add_ships() {
	// Create a new ship mesh
	mesh_drawable shipbody;
	shipbody.initialize(mesh_load_file_obj("assets/shipbody.obj"), "ship");
	shipbody.transform.scaling = 0.0005f;
	shipbody.shading.phong.specular = 0;
	shipbody.shading.phong.diffuse = 0.3f;
	shipbody.shading.phong.ambient = 0.7f;
	shipbody.shading.color = { 0.73f, 0.55f, 0.39f }; // Wood color

	mesh_drawable sail;
	sail.initialize(mesh_load_file_obj("assets/sail.obj"), "sail");
	sail.transform.scaling = 0.0005f;
	sail.shading.phong.specular = 0;
	sail.shading.phong.diffuse = 0.3f;
	sail.shading.phong.ambient = 0.7f;
	sail.shading.color = { 0.996f, 0.996f, 0.886f }; // Wood color

	ship_drawable.add(shipbody);
	ship_drawable.add(sail);

	// Generates random boats on the sea
	for (int i = 0; i < SHIP_CNT; i++) {
		ships.push_back({});
		
		// We add ships wherever there's a ocean biotope
		while (true) {
			int idx = rand() % N;
			if (biotopes[idx] == Biotope::Ocean) {
				vec2 dir = { 0.0f, 1.0f };
				ships[i] = { centers[idx], dir, false, idx };
				break;
			}
		}
	}
}

/// Adds a cloud to the scene.
void scene_structure::add_cloud() {
	// Create a cloud with n*n particles
    normal_distribution<double> distEl4(max_height,0.04f);
	for (int k = 0; k < 3; k++) {
		for (int i = 0; i < PARTICLE_CNT; i++) {
			for (int j = 0; j < PARTICLE_CNT; j++) {
				particles.push_back({(float) i/PARTICLE_CNT, (float) j/PARTICLE_CNT, distEl4(randeng) + k * PARTICLE_SIZE });
			}
		}
	}
	cloud.initialize(mesh_primitive_quadrangle({ -PARTICLE_SIZE,0,0 }, { PARTICLE_SIZE,0,0 }, { PARTICLE_SIZE,0,PARTICLE_SIZE }, { -PARTICLE_SIZE,0,PARTICLE_SIZE }));
	cloud.texture = opengl_load_texture_image("assets/smoke.png");
	cloud.shading.alpha = 0.8f;
	cloud.transform.translation.z = 1.0f;
	cloud.shading.phong.specular = 0.0f;
	cloud.shading.phong.diffuse = 0.0f;
	cloud.shading.phong.ambient = 1.0f;
}

/// Adds snowflakes to the scene.
void scene_structure::add_snowflakes() {
	snowheight = normal_distribution<double>(max_height + 0.7f,0.3f);
	for (int idx = 0; idx < N; idx++) {
		if (biotopes[idx] == Biotope::Snow) {
			vec3 pos = {centers[idx].x, centers[idx].y, snowheight(randeng)};
			snowflakes.push_back({ pos, 0, idx, idx });
		}
	}

	// Create a snowflake mesh
	snowflake.initialize(mesh_primitive_quadrangle({ -PARTICLE_SIZE,0,0 }, { PARTICLE_SIZE,0,0 }, { PARTICLE_SIZE,0,PARTICLE_SIZE }, { -PARTICLE_SIZE,0,PARTICLE_SIZE }));
	snowflake.texture = opengl_load_texture_image("assets/snowflake.png");
	snowflake.shading.phong.specular = 0.0f;
	snowflake.shading.phong.diffuse = 0.0f;
	snowflake.shading.phong.ambient = 1.0f;
}

/// Adds birds/boids to the scene.
void scene_structure::add_birds() {
	// draw birds
	mesh bird_mesh = mesh_load_file_obj("assets/camel.obj");

	for (int i = 0; i < n_birds; i++) {
		birds.push_back({});

		hierarchy_mesh_drawable bird_drawable;

		mesh_drawable body;
		body.initialize(mesh_primitive_ellipsoid({ 1.0f, 1.0f, 1.5f }), "body");
		body.texture = opengl_load_texture_image("assets/body.png"); // associate a texture-image to each element

		mesh_drawable head;
		head.initialize(mesh_primitive_sphere(0.9f), "head");
		head.texture = opengl_load_texture_image("assets/head.png");

		mesh_drawable wing_L;
		wing_L.initialize(mesh_load_file_obj("assets/wing.obj"), "wing_L");
		wing_L.texture = opengl_load_texture_image("assets/wing.png");

		mesh_drawable wing_R;
		wing_R.initialize(mesh_load_file_obj("assets/wing.obj"), "wing_R");
		wing_R.texture = opengl_load_texture_image("assets/wing.png");
		wing_R.transform.rotation = rotation_transform::from_axis_angle({ 0,0,1 }, Pi);

		mesh_drawable eye_L;
		eye_L.initialize(mesh_primitive_sphere(0.2f), "eye_R");
		eye_L.shading.color = { 0.2f,0.2f,0.2f };

		mesh_drawable eye_R;
		eye_R.initialize(mesh_primitive_sphere(0.2f), "eye_L");
		eye_R.shading.color = { 0.2f,0.2f,0.2f };

		bird_drawable.add(body);
		bird_drawable.add(head, "body", { 0,0,1.8 });
		bird_drawable.add(eye_R, "head", { 0.4f,0.3f,0.6f });
		bird_drawable.add(eye_L, "head", { -0.4f,0.3f,0.6f });
		bird_drawable.add(wing_L, "body", { 0,1.03f,0 });
		bird_drawable.add(wing_R, "body", { 0,1.03f,0 });

		vec3 initial_pos = { rand() % 10, rand() % 10, 8 + rand()%2 };
		vec3 initial_speed = vec3( rand_interval(1.0, 10.0), rand_interval(1.0, 10.0), 0 );
		birds[i].pos = initial_pos;
		birds[i].speed = initial_speed;
		birds[i].bird_drawable = bird_drawable;
	}
}


void scene_structure::update_birds_positions(float dt, int n_birds, vector<Bird> &birds) {
	update_birds_speeds(dt, n_birds, birds);
	for (int i = 0; i < n_birds; i++) {
		birds[i].set_position(birds[i].pos + birds[i].speed * dt);
	}
}
void scene_structure::update_birds_speeds(float dt, int n_birds, vector<Bird> &birds) {
	for (int i = 0; i < n_birds; i++) {
		birds[i].fly_towards_others(dt, n_birds, birds);
		birds[i].avoid_collision(dt, n_birds, birds);
		birds[i].adapt_speed_to_others(dt, n_birds, birds);
		birds[i].keep_within_boundaries(dt);
		birds[i].limit_speed(dt);
	}
}
void Bird:: fly_towards_others(float dt, int n_birds, vector<Bird> &birds) {
	float radius_searched = 3.0f;
	vec3 center = { 0,0,0 };
	int n_neigbours = 0;
	for (int i = 0; i < n_birds; i++) {
		if (norm(this->pos - birds[i].pos) < pow(radius_searched, 2)) {
			center += (this->pos - birds[i].pos) * center_follow_factor;
			n_neigbours++;
		}
	}
	if (n_neigbours) this->speed += (center - this->pos) * dt / n_neigbours;
}
void Bird::avoid_collision(float dt, int n_birds, vector<Bird> &birds) {
	float max_distance = 0.8f;
	for (int i = 0; i < n_birds; i++) {
		float dist = norm(this->pos - birds[i].pos);
		if (&birds[i] != this &&  dist < pow(max_distance, 2)) {
			this->speed += (this->pos - birds[i].pos) * dt * avoiding_factor * (pow(dist,2) - 20*dist + 21);
		}
	}
}
void Bird::keep_within_boundaries(float dt) {
	if (pos.x < 0.5) {
		speed.x += border_avoiding_factor;
	}
	if (pos.x > TERRAIN_SIZE) {
		speed.x -= border_avoiding_factor;
	}
	if (pos.y < 0.5) {
		speed.y += border_avoiding_factor;
	}
	if (pos.y > TERRAIN_SIZE) {
		speed.y -= border_avoiding_factor;
	}
	if (pos.z < 5.0) {
		speed.z += border_avoiding_factor;
	}
	if (pos.z > 9.5) {
		speed.z -= border_avoiding_factor;
	}
}
void Bird::adapt_speed_to_others(float dt, int n_birds, vector<Bird>& birds) {
	float radius_searched = 1.0f;
	vec3 neighbours_mean_speed = { 0,0,0 };
	int n_neigbours = 0;
	for (int i = 0; i < n_birds; i++) {
		if (&birds[i] != this && norm(this->pos - birds[i].pos) < pow(radius_searched, 2)) {
			neighbours_mean_speed += birds[i].speed;
			n_neigbours++;
		}
	}
	if (n_neigbours) this->speed += neighbours_mean_speed * adapt_speed_factor / n_neigbours;
}
void Bird::limit_speed(float dt) {
	if (speed.z > 0.5 || speed.z < -0.5) {
		speed.z /= 3.0;
	}
	if (norm(this->speed) > max_speed) {
		this->speed *= max_speed / norm(this->speed);
	}
	else if(norm(this->speed) < min_speed) {
		this->speed *= min_speed / norm(this->speed);
	}
}

void Bird::set_position(vec3 position) {
	this->pos = position;
}

void Bird::set_speed(vec3 speed) {
	this->speed = speed;
}