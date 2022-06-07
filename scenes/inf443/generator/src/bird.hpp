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
#include "wind.hpp"

using namespace cgp;
using namespace std::chrono;
using namespace std;

//A bird
struct Bird {
	hierarchy_mesh_drawable bird_drawable;
	vec3 pos;
	vec3 speed;
	float center_follow_factor = 0.5;
	float adapt_speed_factor = 0.05;
	float avoiding_factor = 3.0;
	float border_avoiding_factor = 0.075;
	float min_speed = 0.5;
	float max_speed = 5.0;
	void fly_towards_others(float dt, int n_birds, vector<Bird> &birds);
	void avoid_collision(float dt, int n_birds, vector<Bird> &birds);
	void keep_within_boundaries(float dt);
	void adapt_speed_to_others(float dt, int n_birds, vector<Bird>& birds);
	void limit_speed(float dt);
	void set_speed(vec3 speed);
	void set_position(vec3 position);
};