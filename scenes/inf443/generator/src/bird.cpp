#include "bird.hpp"

void Bird::fly_towards_others(float dt, int n_birds, vector<Bird> &birds) {
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