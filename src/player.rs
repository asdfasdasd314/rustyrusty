use crate::physics::*;
use raylib::prelude::*;

pub struct Player {
    pub camera: Camera3D,
    pub camera_sensitivity: f32,
    pub movement_speed: f32,
    pub pitch: f32,
    pub yaw: f32,

    pub rigid_body: PhysicsBody,
}

impl Player {
    pub fn new(movement_speed: f32, camera_sensitivity: f32, rigid_body: PhysicsBody) -> Self {
        Player {
            camera: Camera3D::perspective(
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(1.0, 1.0, 1.0),
                Vector3::new(0.0, 1.0, 0.0),
                60.0,
            ),
            movement_speed,
            camera_sensitivity,
            pitch: 0.0,
            yaw: 89.0,
            rigid_body,
        }
    }

    // This function should be called every frame to update the player position
    pub fn update(&mut self, rl: &RaylibHandle) {
        // Update camera position based on input
        if rl.is_key_down(KeyboardKey::KEY_W) {
            self.camera.position.z += self.yaw.to_radians().sin() * self.movement_speed;
            self.camera.position.x += self.yaw.to_radians().cos() * self.movement_speed;
        } else if rl.is_key_down(KeyboardKey::KEY_S) {
            self.camera.position.z -= self.yaw.to_radians().sin() * self.movement_speed;
            self.camera.position.x -= self.yaw.to_radians().cos() * self.movement_speed;
        }
        if rl.is_key_down(KeyboardKey::KEY_A) {
            self.camera.position.z -= self.yaw.to_radians().cos() * self.movement_speed;
            self.camera.position.x += self.yaw.to_radians().sin() * self.movement_speed;
        } else if rl.is_key_down(KeyboardKey::KEY_D) {
            self.camera.position.z += self.yaw.to_radians().cos() * self.movement_speed;
            self.camera.position.x -= self.yaw.to_radians().sin() * self.movement_speed;
        }

        // Also update target vector

        // Get the mouse movement
        let mouse_delta = rl.get_mouse_delta();

        // Update pitch/yaw based on input (replace this with your logic)
        self.yaw += mouse_delta.x * self.camera_sensitivity; // Horizontal rotation
        self.pitch -= mouse_delta.y * self.camera_sensitivity; // Vertical rotation (inverted Y-axis)

        self.pitch = self.pitch.clamp(-89.0, 89.0);

        // Convert angles to direction vector
        let direction = Vector3::new(
            self.yaw.to_radians().cos() * self.pitch.to_radians().cos(),
            self.pitch.to_radians().sin(),
            self.yaw.to_radians().sin() * self.pitch.to_radians().cos(),
        );

        // Update camera target based on direction
        self.camera.target = self.camera.position + direction;
    }
}
