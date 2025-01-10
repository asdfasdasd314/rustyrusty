use raylib::prelude::*;

use crate::game::*;
use crate::physics::*;
use crate::player::*;
use crate::render::*;

mod debug;
mod game;
mod heap;
mod math_util;
mod physics;
mod hashable_vector;
mod player;
mod render;
mod round;
mod world_gen;

fn main() {
    let player_rect = RectangularPrism::new(Vector3::new(0.0, 0.0, 0.0), 1.0, 1.0, 1.0);
    let rb = DynamicBody::new(player_rect.root, Box::new(player_rect));
    let player = Player::new(rb.absolute_position, 10.0, 0.1, CameraType::ThirdPerson(5.0), rb);
    let mut game = Game::new(player, vec![], vec![], (640, 480), "Game".to_string());

    game.game_loop();
}
