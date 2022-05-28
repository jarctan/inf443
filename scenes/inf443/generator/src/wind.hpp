#pragma once

#include "cgp/cgp.hpp"
#include "parameters.hpp"

using namespace cgp;
using namespace std::chrono;
using namespace std;

/// A windsource.
class Windsource {
public:
	virtual ~Windsource() {

	}
	/// Computes the value of the wind at the given position.
	virtual vec2 field(vec2 pos) = 0;
	/// Performs a step in the animation of the source
	/// (it often means moving the source).
	virtual void step(float dt) = 0;
};

/// A gust of wind.
class Gust: public Windsource {
public:
	Gust(vec2 _center, vec2 _direction, float _strength, float _bandwidth)
		: strength(_strength), bandwidth(_bandwidth), direction(_direction), center(_center) { }

	vec2 field(vec2 pos) override {
		float proj = sum((pos - center) * normalize(direction));
		if (proj > - bandwidth / 4 && proj < bandwidth / 4)  {
			return strength * normalize(direction);
		} else if (proj <= - bandwidth / 4) {
			float radius = bandwidth / 4;
			return strength * exp(-0.5f * pow(proj + bandwidth / 4, 2.0f)/pow(radius, 2.0f)) * normalize(direction);
		} else {
			float radius = bandwidth / 4;
			return strength * exp(-0.5f * pow(proj - bandwidth / 4, 2.0f)/pow(radius, 2.0f)) * normalize(direction);
		} 
	}
	void step(float dt) override {
		center += direction * dt;
		if (center.x >= (float) TERRAIN_SIZE || center.x <= 0.0f || center.y >= (float) TERRAIN_SIZE || center.y <= 0.0f) {
			direction = - direction;
		}
	}
private:
	float strength;
	float bandwidth;
	vec2 direction;
	vec2 center;
};

/// A cyclone.
class Cyclone: public Windsource {
public:
	Cyclone(vec2 _center, float _strength, float _radius, vec2 _direction, bool _clockwise = true)
		: center(_center), strength(_strength), radius(_radius), clockwise(_clockwise), direction(_direction) { }

	vec2 field(vec2 pos) override {
		vec2 v = pos - center;
		vec2 t = { - v.y, v.x };
		return (clockwise ? 1 : -1) * strength * exp(-0.5f * pow(norm(v), 2.0f)/pow(radius, 2.0f)) * normalize(t);
	}
	
	void step(float dt) override {
		center += direction * dt;
		if (center.x >= (float) TERRAIN_SIZE || center.x <= 0.0f || center.y >= (float) TERRAIN_SIZE || center.y <= 0.0f) {
			direction = - direction;
		}
	}
private:
	vec2 center;
	float strength;
	float radius;
	bool clockwise;
	vec2 direction;
};

/// A headwind to avoid collisions.
class Revultion: public Windsource {
public:
	Revultion(vec2 _center, float _radius, float _mean, float _strength)
		: center(_center), radius(_radius), mean(_mean), strength(_strength) { }

	vec2 field(vec2 pos) override {
		const vec2 v = pos - center;
		const vec2 none = { 0.0f, 0.0f};
		if (norm(v) > 0.001f)
			return strength * exp(-0.5f * pow(norm(v) - mean, 2.0f)/pow(radius, 2.0f)) * normalize(v);
		else
			return none;
	}
	
	void step(float) override {
		
	}
private:
	vec2 center;
	float radius;
	float mean;
	float strength;
};