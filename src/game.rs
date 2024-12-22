use raylib::prelude::*;

use crate::physics::*;
use crate::player::*;

// This stores the game state and effectively runs everything that has to do with overall game logic
pub struct Game {
    player: Player,

    physical_objects: Vec<RigidBody>,

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
        physical_objects: Vec<RigidBody>,
        window_size: (i32, i32),
        window_title: String,
    ) -> Game {
        let (rl, thread) = raylib::init()
            .size(window_size.0, window_size.1)
            .title(&window_title)
            .build();
        let mut game = Game {
            player,
            physical_objects,
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

            // Begin rendering
            let mut draw_handle = self.raylib_handle.begin_drawing(&self.raylib_thread);

            draw_handle.clear_background(Color::RAYWHITE);

            // This covers everything that is rendered in three dimensions
            {
                let mut draw_handle_3d = draw_handle.begin_mode3D(self.player.camera);

                for object in &self.physical_objects {
                    object.mesh.render(&mut draw_handle_3d);
                }
            }

            draw_handle.draw_fps(10, 10);

            let mut all_objects: Vec<&RigidBody> = Vec::new();

            for physical_object in self.physical_objects.iter() {
                all_objects.push(physical_object);
            }

            all_objects.push(&self.player.rigid_body);

            let colliding_objects = find_colliding_objects(&all_objects);
            if colliding_objects.len() > 0 {
                println!("Colliding");
            }
            else {
                println!("Not colliding");
            }
        }
    }

    pub fn add_physical_object(&mut self, new_object: RigidBody) {
        self.physical_objects.push(new_object);
    }
}

fn find_colliding_objects(objects: &[&RigidBody]) -> Vec<(usize, usize)> {
    let mut colliding_objects: Vec<(usize, usize)> = Vec::new();
    for i in 0..objects.len() - 1 {
        for j in i + 1..objects.len() {
            if rigid_bodies_collide(&objects[i], &objects[j]) {
                colliding_objects.push((i, j));
            }
        }
    }

    colliding_objects
}
