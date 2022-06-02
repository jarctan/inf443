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
#include "environment_camera_head.hpp"
#include "parameters.hpp"
#include "wind.hpp"

using namespace cgp;
using namespace std::chrono;
using namespace std;

/// The element of the GUI that are not already stored in other structures.
struct gui_parameters {
	bool display_frame = true;
	bool display_wireframe = false;
};

/// A biotope.
enum class Biotope {
	Land,
	Ocean,
	Shore,
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

/// A neighbor in the Voronoi diagram.
/// It is characterized by the neighboring polygon,
/// and the vertices of the edge between the two
/// polygons.
struct Neighbor {
  int polygon;
  int vertA;
  int vertB;
};

/// An adjacent vertex in the corresponding Delauney triangulation.
/// It is characterized by the adjacent vertex,
/// and the two polygons on either side of the border between these
/// two vertices.
struct Adjacent {
  int vertex;
  int polyA;
  int polyB;
};

/// A ship.
struct Ship {
	vec2 pos;
	vec2 dir;
	bool outofscope;
	int polygon;
};

/// A snowflake.
struct Snowflake {
	vec3 pos;
	float speed_z;
	int initial_polygon;
	int polygon;
};

/// The structure of the custom scene.
struct scene_structure {
public:
	~scene_structure();
	/// The standard global frame.
	mesh_drawable global_frame;
	/// Standard environment controler.
	scene_environment_camera_head environment;
	/// Storage for inputs status (mouse, keyboard, window dimension).
	inputs_interaction_parameters inputs;

	/// The frame display to be called within the animation loop.
	void display();
	/// Display semi-transparent meshes.
	void display_semi_transparent();
	/// Standard initialization to be called before the animation loop.
	void initialize();
	/// The display of the GUI, also called within the animation loop. 
	void display_gui();

	//variables and methods for the player movement
	float speed = 100;
	float initial_camera_pitch = cgp::Pi / 2.0f;
	float initial_camera_yaw = 0;
	float camera_pitch = 0;
	float camera_yaw = 0;
	float mouseSpeed = 1;
	bool cameraCanMove = true;
	void handleKeyPress(GLFWwindow* window, int key, int action);
	void handleMouseMove(GLFWwindow* window);

private:
	/// Standard GUI element storage.
	gui_parameters gui;

	/// Drawable structures to display the Voronoi diagram.
	mesh_drawable terrain;

	/// Ships to display on the sea: the mesh and their positions.
	mesh_drawable ship_drawable;
	vector<Ship> ships;

	// Cloud particles
	vector<vec3> particles;
	mesh_drawable cloud;

	// Snow
    normal_distribution<double> snowheight;
	vector<Snowflake> snowflakes;
	mesh_drawable snowflake;

	/// Drawable bird
	std::vector<mesh_drawable> birds;
	int n_birds = 10;

	/// Skybox.
	skybox_drawable skybox;

	/// The maximum height of the mountains.
	/// This is being used by the clouds to avoid collision.
	float max_height;

	// Structures representing the Voronoi diagram in memory
	int N; // The number of clusters
	vector<vec3> centers; // Their centers
    vector<vector<int>> mycorners; // The set of the corners of a polygon
    vector<vector<Neighbor>> neighbors; // The set of neighbors of a polygon
    vector<Biotope> biotopes; // The biotope of the polygon
    vector<float> waterdists; // The distance to the water
    vector<Windsource*> windsources; // Wind sources

	int N_corners; // The number of corners
    vector<vec3> corners; // The corners
	vector<vector<int>> touches; // The list of polygons a corner  touches
    vector<vector<Adjacent>> adjacents; // The adjacents corners of a corner
    default_random_engine randeng;

	/// Timer used for the animation.
	timer_basic timer;

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
	/// Adds wind primitives.
	void add_wind();
	/// Adds ships.
	void add_ships();
	/// Adds snowflakes to the scene.
	void add_snowflakes();
	/// Adds cloud particles to the scene.
	void add_cloud();
	/// Adds cloud particles to the scene.
	void add_birds();
};