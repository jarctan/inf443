#pragma once

#include "cgp/cgp.hpp"
#include "environment_camera_head.hpp"

/// The element of the GUI that are not already stored in other structures.
struct gui_parameters {
	bool display_frame = true;
	bool display_wireframe = false;
};


/// The structure of the custom scene.
struct scene_structure {
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //

	/// The standard global frame.
	cgp::mesh_drawable global_frame;
	/// Standard environment controler.
	scene_environment_camera_head environment;
	/// Storage for inputs status (mouse, keyboard, window dimension).
	cgp::inputs_interaction_parameters inputs;

	/// Standard GUI element storage.
	gui_parameters gui;

	void draw_segment(cgp::vec3 const& a, cgp::vec3 const& b);

	// Drawable structures to display the Voronoi diagram
	cgp::mesh_drawable particle_sphere;
	cgp::segments_drawable segment;

	// Structures representing the Voronoi diagram in memory
	// TODO: improve it
	int N; // The number of clusters
	std::vector<cgp::vec3> centers; // Their centers
    std::vector<std::tuple<cgp::vec3,cgp::vec3>> edges; // The edges

	/// Timer used for the animation.
	cgp::timer_basic timer;

	// ****************************** //
	// Functions
	// ****************************** //

	/// Standard initialization to be called before the animation loop.
	void initialize();
	/// The frame display to be called within the animation loop.
	void display();
	/// The display of the GUI, also called within the animation loop. 
	void display_gui();
};

void handleKeyPress();
void handleMouseMove();



