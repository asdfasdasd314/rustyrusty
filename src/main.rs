use raylib::prelude::*;

use crate::game::*;
use crate::physics::*;
use crate::player::*;
use crate::render::*;

mod game;
mod heap;
pub mod math_util;
mod physics;
mod player;
mod render;
mod world_gen;

fn main() {
    let player_rect = RectangularPrism::new(Vector3::new(0.0, 0.0, 0.0), 1.0, 1.0, 1.0);
    //let rect1 = RectangularPrism::new(Vector3::new(2.0, 3.0, 0.0), 10.0, 20.0, 2.0);
    let rect2 = RectangularPrism::new(Vector3::new(5.0, 3.0, 1.0), 12.0, 5.0, 3.0);
    //let solid_body1 = RigidBody::new(Box::new(rect1));
    let solid_body2 = RigidBody::new(rect2.root, Box::new(rect2));
    let rb = RigidBody::new(player_rect.root, Box::new(player_rect));
    let player = Player::new(rb.absolute_position, 10.0, 0.1, rb);
    let mut game = Game::new(player, vec![solid_body2], (640, 480), "Game".to_string());

    game.game_loop();
}
