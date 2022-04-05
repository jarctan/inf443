#pragma once


#include "cgp/cgp.hpp"

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
	cgp::scene_environment_basic_camera_spherical_coords environment;
	/// Storage for inputs status (mouse, keyboard, window dimension).
	cgp::inputs_interaction_parameters inputs;

	/// Standard GUI element storage.
	gui_parameters gui;


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





