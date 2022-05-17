#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <queue>
#include <functional>
#include <random>
#include "cgp/cgp.hpp"

/// The element of the GUI that are not already stored in other structures.
struct gui_parameters {
	bool display_frame = true;
	bool display_wireframe = false;
};


enum class Biotope {
	Land,
	Ocean,
	Lake,
	Snow,
	Tundra,
	Bare,
	Scorched,
	Taiga,
	Shrubland,
	TempDesert,
	TempRainForest,
	TempDeciduousForest,
	Grassland,
	TropicalRainForest,
	TropicalSeasonalForest,
	SubtropicalDesert,
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
	cgp::mesh_drawable sea;
	cgp::mesh_drawable terrain;

	// Structures representing the Voronoi diagram in memory
	// TODO: improve it
	int N; // The number of clusters
	std::vector<cgp::vec3> centers; // Their centers
    std::vector<std::vector<int>> mycorners; // The set of the corners of a polygon
    std::vector<std::vector<std::tuple<int,int,int>>> neighbors; // The neighbors of a polygon
    std::vector<Biotope> biotopes; // The biotope of the cluster
    std::vector<float> waterdists; // The distance to the water

	int N_corners; // The number of corners
    std::vector<cgp::vec3> corners; // The corners
	std::vector<std::vector<int>> touches; // The list of polygons a corner  touches
    std::vector<std::vector<std::tuple<int,int,int>>> adjacents; // The adjacents corners of a corner
    std::default_random_engine randeng;

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
	/// Compares the heights of the two points. Used exclusively for Dijkstra.
	static bool CompareHeights(const std::pair<int,int> &a, const std::pair<int,int> &b);
	/// Creates the terrain.
	void create_terrain();
	/// Computes the heights of the Voronoi polygons.
	void compute_heights();
	/// Computes the smallest distance to the water.
	void compute_waterdists();
	/// Creates a new Voronoi structure with `n` clusters
	/// and updates the class fields accordingly.
	void create_voronoi(int n);
	/// Creates an ocean border around the island.
	void create_ocean_border();
	/// Places the center of the polygons at the mean of the elevation its corners.
	void smooth_centers();
	/// Applies Laplacian smoothing to both the centers and the corners of each polygon.
	void laplacian_smoothing();
	/// Adds remaining biotopes.
	void add_biotopes();
};





