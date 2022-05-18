#include "scene.hpp"
#define JC_VORONOI_IMPLEMENTATION
#include "jc_voronoi.h"

using namespace cgp;

/// Place the points at the centers (centroid) of their respective clusters
/// in the Voronoi diagram.
void relax_points(const jcv_diagram* diagram, jcv_point* points) {
    const jcv_site* sites = jcv_diagram_get_sites(diagram);
    for( int i = 0; i < diagram->numsites; ++i )
    {
        const jcv_site* site = &sites[i];
        jcv_point sum = site->p;
        int count = 1;

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
		else if (key == GLFW_KEY_LEFT_SHIFT && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			displacement = speed * dt * vec3(0, 0, -1);
		}

		//std::cout << camera.position() << "\n";
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

		camera_yaw += (screenWidth / 2.0f - xpos) * dt / 50.0f;
		camera_pitch += (screenHeight / 2.0f - ypos) * dt / 50.0f;

		//std::cout << "yaw = " << camera_yaw << " pitch = " << camera_pitch << "\n";

		rotation_transform r_pitch = rotation_transform::from_axis_angle(vec3(1, 0, 0), camera_pitch);
		rotation_transform r_yaw = rotation_transform::from_axis_angle(vec3(0, 1, 0), camera_yaw);

		camera.orientation_camera = initialOrientation * r_yaw * r_pitch;
	}
}

/// This function is called only once at the beginning of the program
/// and initializes the meshes, diagrams and other structures.
void scene_structure::initialize() {
	// Set the behavior of the camera and its initial position
	environment.camera.position_camera = { 5.0f, 5.0f, 10.0f };
	environment.camera.manipulator_rotate_roll_pitch_yaw(0, camera_pitch, camera_yaw); //initial rotation value

	// Number of clusters
	N = 1000;

	jcv_diagram diagram;
    memset(&diagram, 0, sizeof(jcv_diagram));
	jcv_point* points = (jcv_point*) malloc(N*sizeof(jcv_point));
	for (int i = 0; i < N; i++) {
		jcv_point pt;
		pt.x = rand() % 300;
		pt.y = rand() % 300;
		points[i] = pt;
	}

	// Run n iterations of Voronoi (find Voronoi clusters, then place points at centroids, repeat)
	for (int i = 0; i < 2; i++) {
		jcv_diagram_generate(N, points, 0, 0, &diagram );
		relax_points(&diagram, points);
	}

	segments_drawable::default_shader = curve_drawable::default_shader;
	segment.initialize({ {0,0,0},{1,0,0} });

	// Store the centers
    for (int i = 0; i < N; i++) {
		centers.push_back({points[i].x, points[i].y, 0});
	}

	// Store the edges
    const jcv_edge* edge = jcv_diagram_get_edges( &diagram );
    while (edge) {
		vec3 p1 = {edge->pos[0].x, edge->pos[0].y, 0 };
		vec3 p2 = {edge->pos[1].x, edge->pos[1].y, 0 };
        edges.push_back(std::make_tuple(p1, p2));
        edge = jcv_diagram_get_next_edge(edge);
    }

	// Create a mesh for the centers
	particle_sphere.initialize(mesh_primitive_sphere(1.0f));
	particle_sphere.shading.color = { 1,0,0 };

	// Always free up memory at the end
    jcv_diagram_free( &diagram );
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

	for(std::tuple<vec3, vec3>& edge: edges)
		draw_segment(std::get<0>(edge), std::get<1>(edge));

	for (vec3& center : centers) {
		//std::cout << "New center(" << center.x << "; " << center.y << "; " << center.z << ")" << std::endl;
		particle_sphere.transform.translation = center;
		draw(particle_sphere, environment);
	}

	if (gui.display_wireframe);
}

void scene_structure::display_gui() {
	ImGui::Checkbox("Frame", &gui.display_frame);
	ImGui::Checkbox("Wireframe", &gui.display_wireframe);
}
