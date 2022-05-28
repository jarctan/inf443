#pragma once

const float MOIS_STD = 0.2f;
const float EL_STD = 0.2f;
/// Number of relaxations to perform on the Voronoi diagram.
const int RELAX_CNT = 2;
/// (Estimate) number of polygons to use to generate the diagram.
const int POLYGON_CNT = 20000;
/// SIze of the generated terrain (height and width)
const float TERRAIN_SIZE = 20.0f;
/// Square root of the number of particle par cloud layer.
/// IMPORTANT: Change this value if your fps rate is too low.
const int PARTICLE_CNT = 20;
/// Number of ships.
const int SHIP_CNT = 100;