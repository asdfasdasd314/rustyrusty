use crate::physics::*;
use crate::player::*;
use raylib::prelude::*;

mod game;
mod math_util;
mod physics;
mod player;
mod world_gen;

const MAX_COLUMNS: usize = 20;

fn main() {
    let shape = RectangularPrism::new(Vector3::new(10.0, 0.0, 0.0), 10.0, 10.0, 10.0);
    let rigid_body = PhysicsBody::new(Box::new(shape));
    let mut player = Player::new(0.01, 0.1, rigid_body);
    let (mut rl, thread) = raylib::init().size(640, 480).title("Physics!!!").build();

    rl.hide_cursor();
    rl.disable_cursor();

    // Generates some random columns
    let mut heights = [0.032; MAX_COLUMNS];
    let mut positions = [Vector3::zero(); MAX_COLUMNS];
    let mut colors = [Color::default(); MAX_COLUMNS];

    for i in 0..MAX_COLUMNS {
        heights[i] = rl.get_random_value::<i32>(1..12) as f32;
        positions[i] = rvec3(
            rl.get_random_value::<i32>(-15..15),
            heights[i] / 2.0,
            rl.get_random_value::<i32>(-15..15),
        );
        colors[i] = Color::new(
            rl.get_random_value::<i32>(20..255) as u8,
            rl.get_random_value::<i32>(10..55) as u8,
            30,
            255,
        );
    }

    let mut cursor_shown = false;

    while !rl.window_should_close() {
        if rl.is_key_pressed(KeyboardKey::KEY_ENTER) {
            rl.toggle_fullscreen();
        }
        if rl.is_key_pressed(KeyboardKey::KEY_TAB) {
            if cursor_shown {
                rl.hide_cursor();
                rl.disable_cursor();
            } else {
                rl.show_cursor();
            }
            cursor_shown = !cursor_shown;
        }

        // Update
        //----------------------------------------------------------------------------------
        player.update(&rl);
        //----------------------------------------------------------------------------------

        // Draw
        //----------------------------------------------------------------------------------
        let mut d = rl.begin_drawing(&thread);

        d.clear_background(Color::RAYWHITE);

        {
            let mut d = d.begin_mode3D(player.camera);

            d.draw_plane(rvec3(0.0, 0.0, 0.0), rvec2(32.0, 32.0), Color::LIGHTGRAY); // Draw ground
            d.draw_cube(rvec3(16.0, 2.5, 0.0), 1.0, 5.0, 32.0, Color::BLUE); // Draw a blue wall
            d.draw_cube(rvec3(16.0, 2.5, 0.0), 1.0, 5.0, 32.0, Color::LIME); // Draw a green wall
            d.draw_cube(rvec3(0.0, 2.5, 16.0), 32.0, 5.0, 1.0, Color::GOLD); // Draw a yellow wall

            // Draw some cubes around
            for i in 0..MAX_COLUMNS {
                d.draw_cube(positions[i], 2.0, heights[i], 2.0, colors[i]);
                d.draw_cube_wires(positions[i], 2.0, heights[i], 2.0, Color::MAROON);
            }
        }

        d.draw_rectangle(10, 10, 220, 70, Color::SKYBLUE.fade(0.5));
        d.draw_rectangle_lines(10, 10, 220, 70, Color::BLUE);

        d.draw_text(
            "First person camera default controls:",
            20,
            20,
            10,
            Color::BLACK,
        );
        d.draw_text("- Move with keys: W, A, S, D", 40, 40, 10, Color::DARKGRAY);
        d.draw_text("- Mouse move to look around", 40, 60, 10, Color::DARKGRAY);

        //----------------------------------------------------------------------------------   while !rl.window_should_close() {
    }
}
