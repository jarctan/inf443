#pragma once

/// Standard deviation of the border btw the moisture zones.
const float MOIS_STD = 0.2f;
/// Standard deviation of the border btw the elevation zones.
const float EL_STD = 0.2f;
/// Number of relaxations to perform on the Voronoi diagram.
const int RELAX_CNT = 2;
/// (Estimate) number of polygons to use to generate the diagram.
const int POLYGON_CNT = 15000;
/// SIze of the generated terrain (height and width)
const float TERRAIN_SIZE = 20.0f;
/// Square root of the number of particle par cloud layer.
/// IMPORTANT: Change this value if your fps rate is too low.
const int PARTICLE_CNT = 20;
/// Number of ships.
/// IMPORTANT: Change this value if your fps rate is too low.
const int SHIP_CNT = 40;
/// Size of a particle (whether it be a cloud or a snowflake).
const float PARTICLE_SIZE = 0.05f;