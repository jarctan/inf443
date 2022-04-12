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

/// This function is called only once at the beginning of the program
/// and initializes the meshes, diagrams and other structures.
void scene_structure::initialize() {
	// Set the behavior of the camera and its initial position
	environment.camera.axis = camera_spherical_coordinates_axis::z;
	environment.camera.look_at({ -1.0f, 4.0f, 2.0f } /*eye position*/, {0,0,0} /*target position*/);

	// Create a visual frame representing the coordinate system
	global_frame.initialize(mesh_primitive_frame(), "Frame");
	environment.camera.axis = camera_spherical_coordinates_axis::z;
	environment.camera.look_at({ 2.0f,-2.0f,1.0f }, { 0,0,0 });

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
