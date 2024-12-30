use raylib::prelude::*;

use crate::physics::*;
use crate::player::*;

// This stores the game state and effectively runs everything that has to do with overall game logic
pub struct Game {
    player: Player,

    dynamic_objects: Vec<RigidBody>,
    static_objects: Vec<RigidBody>,

    window_size: (i32, i32),
    window_title: String,
    cursor_shown: bool,
    fullscreen: bool,

    raylib_handle: RaylibHandle,
    raylib_thread: RaylibThread,
}

impl Game {
    // Handles all the setup for a game (window, player, world, etc.)
    pub fn new(
        player: Player,
        dynamic_objects: Vec<RigidBody>,
        static_objects: Vec<RigidBody>,
        window_size: (i32, i32),
        window_title: String,
    ) -> Game {
        let (rl, thread) = raylib::init()
            .size(window_size.0, window_size.1)
            .title(&window_title)
            .build();
        let mut game = Game {
            player,
            dynamic_objects,
            static_objects,
            window_size,
            window_title,
            cursor_shown: false,
            fullscreen: false,
            raylib_handle: rl,
            raylib_thread: thread,
        };

        game.raylib_handle.hide_cursor();
        game.raylib_handle.disable_cursor();
        return game;
    }

    // Perhaps this could change in the future, but the game loop loops forever and never returns
    pub fn game_loop(&mut self) {
        while !self.raylib_handle.window_should_close() {
            // Calculate delta time
            let delta_time = self.raylib_handle.get_frame_time();
            // When it comes to the game loop there are a few parts: process inputs, create outputs, render

            // Do things about fullscreen and exiting
            if self.raylib_handle.is_key_pressed(KeyboardKey::KEY_ENTER) {
                self.raylib_handle.toggle_fullscreen();
            }
            if self.raylib_handle.is_key_pressed(KeyboardKey::KEY_TAB) {
                if self.cursor_shown {
                    self.raylib_handle.hide_cursor();
                    self.raylib_handle.disable_cursor();
                } else {
                    self.raylib_handle.show_cursor();
                }
                self.cursor_shown = !self.cursor_shown;
            }

            // Do user input and movement
            self.player.update(&self.raylib_handle, delta_time);

            let colliding_objects: Vec<(&RigidBody, &RigidBody)> = self.find_colliding_objects();

            // Begin rendering
            let mut draw_handle = self.raylib_handle.begin_drawing(&self.raylib_thread);

            draw_handle.clear_background(Color::RAYWHITE);

            // This covers everything that is rendered in three dimensions
            {
                let mut draw_handle_3d = draw_handle.begin_mode3D(self.player.camera);

                for object in &self.static_objects {
                    object.mesh.render(&mut draw_handle_3d);
                }

                for object in &self.dynamic_objects {
                    object.mesh.render(&mut draw_handle_3d);
                }

                match self.player.camera_type {
                    CameraType::FirstPerson => {}
                    CameraType::ThirdPerson(_) => {
                        self.player.rigid_body.mesh.render(&mut draw_handle_3d);
                    }
                }
            }

            draw_handle.draw_fps(10, 10);
        }
    }

    pub fn add_dynamic_object(&mut self, new_object: RigidBody) {
        self.dynamic_objects.push(new_object);
    }

    pub fn add_static_object(&mut self, new_object: RigidBody) {
        self.static_objects.push(new_object);
    }

    fn find_colliding_objects(&self) -> Vec<(&RigidBody, &RigidBody)> {
        let mut dynamic_objects: Vec<&RigidBody> = Vec::with_capacity(self.dynamic_objects.len() + 1);
        dynamic_objects.push(&self.player.rigid_body);
        for object in self.dynamic_objects.iter() {
            dynamic_objects.push(object);
        }

        let mut colliding_objects: Vec<(&RigidBody, &RigidBody)> = Vec::new();
        for i in 0..dynamic_objects.len() {
            let dynamic_object1 = dynamic_objects[i];
            let radius1 = dynamic_object1.get_bounding_circle_radius();
            for static_object in self.static_objects.iter() {
                let distance_between_centers = (dynamic_object1.get_center() - static_object.get_center()).length();
                if radius1 + static_object.get_bounding_circle_radius() < distance_between_centers {
                    continue;
                }

                if dynamic_object1.collides_with(static_object) {
                    colliding_objects.push((dynamic_object1, static_object));
                }
            }

            for j in i + 1..dynamic_objects.len() {
                let dynamic_object2 = dynamic_objects[j];
                let distance_between_centers = (dynamic_object1.get_center() - dynamic_object2.get_center()).length();
                if radius1 + dynamic_object2.get_bounding_circle_radius() < distance_between_centers {
                    continue;
                }
                if dynamic_object1.collides_with(dynamic_object2) {
                    colliding_objects.push((dynamic_object1, dynamic_object2));
                }
            }
        }

        colliding_objects
    }
}
