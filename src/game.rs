use crate::physics::*;
use crate::player::*;
use raylib::prelude::*;

// These usizes are the indices in the dynamic and static object vecs in Game struct
enum CollisionObject1 {
    Player,
    Dynamic(usize),
}

enum CollisionObject2 {
    Dynamic(usize),
    Static(usize),
}

// This stores the game state and effectively runs everything that has to do with overall game logic
pub struct Game {
    player: Player,

    dynamic_objects: Vec<DynamicBody>,
    static_objects: Vec<StaticBody>,

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
        dynamic_objects: Vec<DynamicBody>,
        static_objects: Vec<StaticBody>,
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

            let collisions: Vec<(CollisionObject1, CollisionObject2, Vector3)> = self.find_colliding_objects();
            if collisions.len() > 0 {
                self.simulate_collisions(collisions);
            }

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
                        self.player.dynamic_body.mesh.render(&mut draw_handle_3d);
                    }
                }
            }

            draw_handle.draw_fps(10, 10);
        }
    }

    pub fn add_dynamic_object(&mut self, new_object: DynamicBody) {
        self.dynamic_objects.push(new_object);
    }

    pub fn add_static_object(&mut self, new_object: StaticBody) {
        self.static_objects.push(new_object);
    }

    /**
    Returns the first collision object that is said to have "initialized" the collision (either a player or other dynamic body in the world)
    Second entry in the return is the second collision that is "hit" by the first collision object (either a dynamic object or static body, and if it were multiplayer it could be another player object)
    Third entry is the mtv which is applied to the first object to fix intersections between objects
     */
    fn find_colliding_objects(&self) -> Vec<(CollisionObject1, CollisionObject2, Vector3)> {
        let mut collisions: Vec<(CollisionObject1, CollisionObject2, Vector3)> = Vec::new();
        // Check player collisions
        let player_object = &self.player.dynamic_body;
        let player_radius = player_object.get_bounding_circle_radius();
        for i in 0..self.dynamic_objects.len() {
            let dynamic_object = &self.dynamic_objects[i];
            let distance_between_centers = (player_object.get_center() - dynamic_object.get_center()).length();
            if player_radius + dynamic_object.get_bounding_circle_radius() < distance_between_centers {
                continue;
            }
            let collided = player_object.collides_with(Box::new(dynamic_object as &dyn Physical));
            match collided {
                Some(mtv) => {
                    let collision_object1 = CollisionObject1::Player;
                    let collision_object2 = CollisionObject2::Dynamic(i);
                    collisions.push((collision_object1, collision_object2, mtv));
                }
                None => {}
            }
        }
        for i in 0..self.static_objects.len() {
            let static_object = &self.static_objects[i];
            let distance_between_centers = (player_object.get_center() - static_object.get_center()).length();
            if player_radius + static_object.get_bounding_circle_radius() < distance_between_centers {
                continue;
            }
            let collided = player_object.collides_with(Box::new(static_object as &dyn Physical));
            match collided {
                Some(mtv) => {
                    let collision_object1 = CollisionObject1::Player;
                    let collision_object2 = CollisionObject2::Static(i);
                    collisions.push((collision_object1, collision_object2, mtv));
                }
                None => {}
            }
        }

        // Check all other dynamic object collisions
        for i in 0..self.dynamic_objects.len() {
            let dynamic_object1 = &self.dynamic_objects[i];
            let radius1 = dynamic_object1.get_bounding_circle_radius();
            for j in 0..self.static_objects.len() {
                let static_object = &self.static_objects[j];
                let distance_between_centers =
                    (dynamic_object1.get_center() - static_object.get_center()).length();
                if radius1 + static_object.get_bounding_circle_radius() < distance_between_centers {
                    continue;
                }

                let collided = dynamic_object1.collides_with(Box::new(static_object as &dyn Physical));
                match collided {
                    Some(mtv) => {
                        let collision_object1 = CollisionObject1::Dynamic(i);
                        let collision_object2 = CollisionObject2::Static(j);
                        collisions.push((collision_object1, collision_object2, mtv));
                    }
                    None => {}
                }
            }

            for j in i + 1..self.dynamic_objects.len() {
                let dynamic_object2 = &self.dynamic_objects[j];
                let distance_between_centers =
                    (dynamic_object1.get_center() - dynamic_object2.get_center()).length();
                if radius1 + dynamic_object2.get_bounding_circle_radius() < distance_between_centers
                {
                    continue;
                }

                let collided = dynamic_object1.collides_with(Box::new(dynamic_object2 as &dyn Physical));
                match collided {
                    Some(mtv) => {
                        let collision_object1 = CollisionObject1::Dynamic(i);
                        let collision_object2 = CollisionObject2::Dynamic(j);
                        collisions.push((collision_object1, collision_object2, mtv));
                    }
                    None => {}
                }
            }
        }

        collisions
    }

    fn simulate_collisions(&mut self, collisions: Vec<(CollisionObject1, CollisionObject2, Vector3)>) {
        for collision in collisions {
            println!("{:#?}", collision.2);
            match collision.0 {
                CollisionObject1::Player => {
                    self.player.move_by(collision.2);
                }
                CollisionObject1::Dynamic(index) => {
                    self.dynamic_objects[index].move_by(collision.2);
                }
            }
        }
    }
}
