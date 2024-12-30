use raylib::prelude::*;

use crate::game::*;
use crate::physics::*;
use crate::player::*;
use crate::render::*;

mod game;
mod heap;
mod math_util;
mod physics;
mod hashable_vector;
mod player;
mod render;
mod world_gen;

fn main() {
    let mut tons_of_objects: Vec<RigidBody> = Vec::with_capacity(100);
    for i in 100..110 {
        let obj = RectangularPrism::new(Vector3::new(i as f32, i as f32, i as f32), 0.5, 0.5, 0.5);
        let rb = RigidBody::new(obj.root, Box::new(obj));
        tons_of_objects.push(rb);
    }

    let player_rect = RectangularPrism::new(Vector3::new(0.0, 0.0, 0.0), 1.0, 1.0, 1.0);
    //let rect1 = RectangularPrism::new(Vector3::new(2.0, 3.0, 0.0), 10.0, 20.0, 2.0);
    let rect2 = RectangularPrism::new(Vector3::new(5.0, 3.0, 1.0), 12.0, 5.0, 3.0);
    //let solid_body1 = RigidBody::new(Box::new(rect1));
    let solid_body2 = RigidBody::new(rect2.root, Box::new(rect2));
    let rb = RigidBody::new(player_rect.root, Box::new(player_rect));
    //let player = Player::new(rb.absolute_position, 10.0, 0.1, CameraType::ThirdPerson(5.0), rb);
    let player = Player::new(rb.absolute_position, 10.0, 0.1, CameraType::ThirdPerson(5.0), rb);
    let mut game = Game::new(player, tons_of_objects, vec![solid_body2], (640, 480), "Game".to_string());

    game.game_loop();
}
