#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <queue>
#include <functional>
#include "cgp/cgp.hpp"

/// The element of the GUI that are not already stored in other structures.
struct gui_parameters {
	bool display_frame = true;
	bool display_wireframe = false;
};


enum class Biotope {
	Land,
	Ocean
};

/// The structure of the custom scene.
struct scene_structure {
public:
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

	void draw_segment(cgp::vec3 const& a, cgp::vec3 const& b);

	// Drawable structures to display the Voronoi diagram
	cgp::mesh_drawable particle_sphere;
	cgp::segments_drawable segment;

	// Structures representing the Voronoi diagram in memory
	// TODO: improve it
	int N; // The number of clusters
	std::vector<cgp::vec3> centers; // Their centers
    std::vector<cgp::vec3> corners; // The corners
    std::vector<std::vector<std::tuple<int,int,int>>> neighbors; // The neighbors
	std::vector<std::vector<int>> touches; // The list of polygons a cluster touches
    std::vector<Biotope> biotopes; // The biotope of the cluster
    std::vector<float> heights; // The elevation of the cluster

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
private:
	/// Compare the heights of the two points. Used exclusively for Dijkstra
	static bool CompareHeights(const std::pair<int,int> &a, const std::pair<int,int> &b);
};





