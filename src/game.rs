use crate::physics::*;
use crate::player::*;
use raylib::prelude::*;

// This stores the game state and effectively runs everything that has to do with overall game logic
pub struct Game {
    player: Player,
    physical_objects: Vec<SolidBody>,

    window_size: (i32, i32),
    window_title: String,
    cursor_shown: bool,
    fullscreen: bool,

    raylib_handle: RaylibHandle,
    raylib_thread: RaylibThread,
}

impl Game {
    // Handles all the setup for a game (window, player, world, etc.)
    pub fn new(player: Player, physical_objects: Vec<SolidBody>, window_size: (i32, i32), window_title: String) -> Game {
        let (rl, thread) = raylib::init().size(window_size.0, window_size.1).title(&window_title).build();
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
            // When it comes to the game loop there are a few parts: process inputs, create outputs, render

            // Do things about fullscreen and exiting
            if self.raylib_handle.is_key_pressed(KeyboardKey::KEY_ENTER) {
                self.raylib_handle.toggle_fullscreen();
            }
            if self.raylib_handle.is_key_pressed(KeyboardKey::KEY_TAB) {
                if self.cursor_shown {
                    self.raylib_handle.hide_cursor();
                    self.raylib_handle.disable_cursor();
                }
                else {
                    self.raylib_handle.show_cursor();
                }
               self.cursor_shown = !self.cursor_shown;
            }
        

            // Do user input and movement
            self.player.update(&self.raylib_handle);

            // Begin rendering
            let mut d = self.raylib_handle.begin_drawing(&self.raylib_thread);

            d.clear_background(Color::RAYWHITE);

            // This covers everything that is rendered in three dimensions
            {
                let mut d = d.begin_mode3D(self.player.camera);

                for object in &self.physical_objects {
                    object.mesh.render(&mut d);
                }
            }
        }
    }

    pub fn add_physical_object(&mut self, new_object: SolidBody) {
        self.physical_objects.push(new_object);
    }
}
