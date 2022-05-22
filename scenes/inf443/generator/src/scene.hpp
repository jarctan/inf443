#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <queue>
#include <functional>
#include <random>
#include <chrono>
#include "cgp/cgp.hpp"
#include "voronoi/include/Point2.h"
#include "voronoi/include/Vector2.h"
#include "voronoi/include/VoronoiDiagramGenerator.h"


using namespace cgp;
using namespace std::chrono;
using namespace std;

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

struct Neighbor {
  int polygon;
  int vertA;
  int vertB;
};

struct Adjacent {
  int vertex;
  int polyA;
  int polyB;
};

/// The structure of the custom scene.
struct scene_structure {
public:
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //

	/// The standard global frame.
	mesh_drawable global_frame;
	/// Standard environment controler.
	scene_environment_basic_camera_spherical_coords environment;
	/// Storage for inputs status (mouse, keyboard, window dimension).
	inputs_interaction_parameters inputs;

	/// Standard GUI element storage.
	gui_parameters gui;

	/// Drawable structures to display the Voronoi diagram.
	mesh_drawable terrain;

	/// Drawable structures to display the Voronoi diagram.
	mesh_drawable ship;

	/// Skybox.
	skybox_drawable skybox;

	// Structures representing the Voronoi diagram in memory
	int N; // The number of clusters
	vector<vec3> centers; // Their centers
    vector<vector<int>> mycorners; // The set of the corners of a polygon
    vector<vector<Neighbor>> neighbors; // The set of neighbors of a polygon
    vector<Biotope> biotopes; // The biotope of the polygon
    vector<float> waterdists; // The distance to the water
    vector<vector<pair<float,float>>> windfield; // Wind field
    vector<vector<float>> heightfield; // Height field

	int N_corners; // The number of corners
    vector<vec3> corners; // The corners
	vector<vector<int>> touches; // The list of polygons a corner  touches
    vector<vector<Adjacent>> adjacents; // The adjacents corners of a corner
    default_random_engine randeng;

	/// Timer used for the animation.
	timer_basic timer;

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
	/// Adds a wind field
	void add_wind();
};