#include "movement.hpp"

void handleKeyPress() {
	
	inputs_keyboard_parameters const& keyboard = inputs.keyboard;
	camera_head& camera = environment.camera;

	// The camera moves forward all the time
	//   We consider in this example a constant velocity, so the displacement is: velocity * dt * front-camera-vector
	float const dt = timer.update();
	vec3 const forward_displacement = gui.speed * 0.1f * dt * camera.front();
	camera.position_camera += forward_displacement;

	// The camera rotates if we press on the arrow keys
	//  The rotation is only applied to the roll and pitch degrees of freedom.
	float const pitch = 0.5f; // speed of the pitch
	float const roll = 0.7f; // speed of the roll
	if (keyboard.up)
		camera.manipulator_rotate_roll_pitch_yaw(0, -pitch * dt, 0);
	if (keyboard.down)
		camera.manipulator_rotate_roll_pitch_yaw(0, pitch * dt, 0);
	if (keyboard.right)
		camera.manipulator_rotate_roll_pitch_yaw(roll * dt, 0, 0);
	if (keyboard.left)
		camera.manipulator_rotate_roll_pitch_yaw(-roll * dt, 0, 0);
}

void handleMouseMove() {

}