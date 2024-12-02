use raylib::prelude::*;

const MAX_COLUMNS: usize = 20;

fn main() {
    let mut camera: Camera3D = Camera3D::perspective(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0), Vector3::new(0.0, 1.0, 0.0), 60.0);
    let (mut rl, thread) = raylib::init()
        .size(640, 480)
        .title("Hello, World")
        .build();

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
    while !rl.window_should_close() {
            // Update
        //----------------------------------------------------------------------------------
        rl.update_camera(&mut camera, CameraMode::CAMERA_FIRST_PERSON);                  // Update camera
        //----------------------------------------------------------------------------------

        // Draw
        //----------------------------------------------------------------------------------
        let mut d = rl.begin_drawing(&thread);

            d.clear_background(Color::RAYWHITE);

            {
                let mut d = d.begin_mode3D(&camera);
    
                    d.draw_plane(rvec3(0.0, 0.0, 0.0), rvec2(32.0, 32.0), Color::LIGHTGRAY); // Draw ground
                    d.draw_cube(rvec3(16.0, 2.5, 0.0), 1.0, 5.0, 32.0, Color::BLUE);     // Draw a blue wall
                    d.draw_cube(rvec3(16.0, 2.5, 0.0), 1.0, 5.0, 32.0, Color::LIME);      // Draw a green wall
                    d.draw_cube(rvec3(0.0, 2.5, 16.0), 32.0, 5.0, 1.0, Color::GOLD);      // Draw a yellow wall
    
                    // Draw some cubes around
                    for i in 0..MAX_COLUMNS
                    {
                        d.draw_cube(positions[i], 2.0, heights[i], 2.0, colors[i]);
                        d.draw_cube_wires(positions[i], 2.0, heights[i], 2.0, Color::MAROON);
                    }

            }


            d.draw_rectangle( 10, 10, 220, 70, Color::SKYBLUE.fade(0.5));
            d.draw_rectangle_lines( 10, 10, 220, 70, Color::BLUE);

            d.draw_text("First person camera default controls:", 20, 20, 10, Color::BLACK);
            d.draw_text("- Move with keys: W, A, S, D", 40, 40, 10, Color::DARKGRAY);
            d.draw_text("- Mouse move to look around", 40, 60, 10, Color::DARKGRAY);

        //----------------------------------------------------------------------------------   while !rl.window_should_close() {
   }
}
