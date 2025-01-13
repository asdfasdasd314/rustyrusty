use raylib::prelude::*;
use crate::float_precision::*;

use crate::physics::*;

pub enum CameraType {
    FirstPerson,
    ThirdPerson(f64), // Gives the distance of which the third person camera is from the center of the player object
}

pub struct Player {
    // This is an absolute point
    pub absolute_position: Vector3f64,
    pub camera: Camera3D,
    pub camera_type: CameraType,
    pub camera_sensitivity: f64,
    pub movement_speed: f64,
    pub pitch: f64,
    pub yaw: f64,

    pub dynamic_body: DynamicBody,
}

impl Dynamic for Player {
    fn move_by(&mut self, change: Vector3f64) {
        self.absolute_position += change;
        self.dynamic_body.move_by(change);
    }
}

impl Player {
    pub fn new(
        init_position: Vector3f64,
        movement_speed: f64,
        camera_sensitivity: f64,
        camera_type: CameraType,
        dynamic_body: DynamicBody,
    ) -> Self {
        Player {
            absolute_position: init_position,
            camera: Camera3D::perspective(
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(1.0, 1.0, 1.0),
                Vector3::new(0.0, 1.0, 0.0),
                60.0,
            ),
            camera_type,
            movement_speed,
            camera_sensitivity,
            pitch: 0.0,
            yaw: 89.0,
            dynamic_body,
        }
    }

    // This function should be called every frame to update the player position
    pub fn update(&mut self, rl: &RaylibHandle, delta_time: f64) {
        // Update camera position based on input
        let mut position_change: Vector3f64 = Vector3f64::new(0.0, 0.0, 0.0);
        if rl.is_key_down(KeyboardKey::KEY_W) {
            position_change.z += self.yaw.to_radians().sin();
            position_change.x += self.yaw.to_radians().cos();
        } else if rl.is_key_down(KeyboardKey::KEY_S) {
            position_change.z += -self.yaw.to_radians().sin();
            position_change.x += -self.yaw.to_radians().cos();
        }
        if rl.is_key_down(KeyboardKey::KEY_A) {
            position_change.z += -self.yaw.to_radians().cos();
            position_change.x += self.yaw.to_radians().sin();
        } else if rl.is_key_down(KeyboardKey::KEY_D) {
            position_change.z += self.yaw.to_radians().cos();
            position_change.x += -self.yaw.to_radians().sin();
        }
        if rl.is_key_down(KeyboardKey::KEY_SPACE) {
            position_change.y += 1.0;
        } else if rl.is_key_down(KeyboardKey::KEY_LEFT_SHIFT) {
            position_change.y += -1.0;
        }
        position_change *= delta_time * self.movement_speed;

        self.move_by(position_change);

        // Also update target vector

        // Get the mouse movement
        let mouse_delta = rl.get_mouse_delta();

        // Update pitch/yaw based on input (replace this with your logic)
        self.yaw += mouse_delta.x as f64 * self.camera_sensitivity; // Horizontal rotation
        self.pitch -= mouse_delta.y as f64 * self.camera_sensitivity; // Vertical rotation (inverted Y-axis)

        self.pitch = self.pitch.clamp(-89.0, 89.0);

        // Convert angles to direction vector
        let direction = Vector3f64::new(
            self.yaw.to_radians().cos() * self.pitch.to_radians().cos(),
            self.pitch.to_radians().sin(),
            self.yaw.to_radians().sin() * self.pitch.to_radians().cos(),
        );


        match &self.camera_type {
            CameraType::FirstPerson => {
                // Update camera target based on direction
                self.camera.target = self.camera.position + Vector3::from(direction);
            }
            CameraType::ThirdPerson(distance) => {
                // Given the direcction, the camera should be on the opposite side of the sphere from the camera direction
                let target_position = self.dynamic_body.get_center() + direction.normalized() * (-1.0 * distance);
                self.camera.position = self.camera.position.lerp(target_position.into(), 1.0); // Adjust the factor (0.1) for smoothness
                self.camera.target = self.dynamic_body.get_center().into();
            }
        }
    }
}
